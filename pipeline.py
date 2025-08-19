#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, re, json, argparse, sys
import requests
import pandas as pd
from cyvcf2 import VCF
import pysam

# ---------- Helpers ----------

CS_MAP = {
    "benign":"Benign",
    "likely_benign":"Likely benign", "likely benign":"Likely benign",
    "pathogenic":"Pathogenic",
    "likely_pathogenic":"Likely pathogenic", "likely pathogenic":"Likely pathogenic",
    "uncertain_significance":"Uncertain significance", "uncertain significance":"Uncertain significance", "vus":"Uncertain significance",
    "conflicting_interpretations_of_pathogenicity":"Conflicting interpretations",
    "conflicting interpretations of pathogenicity":"Conflicting interpretations",
}

def log(msg): print(msg, flush=True)

def sql_escape(val: str) -> str:
    return str(val).replace("\\", "\\\\").replace("'", "\\'")

def make_synth_id(key_tuple):
    if not isinstance(key_tuple, tuple) or len(key_tuple) != 4: return ""
    contig, pos, ref, alts = key_tuple
    alts_str = ",".join(alts) if isinstance(alts, tuple) else str(alts)
    return f"{contig}:{pos}{ref}>{alts_str}"

def download_if_needed(url, path, skip_download=False):
    if os.path.exists(path):
        log(f"✓ Já existe: {path}")
        return
    if skip_download:
        raise FileNotFoundError(f"Arquivo ausente e SKIP_DOWNLOAD=True: {path}")
    log(f"↓ Baixando: {url}")
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(path, "wb") as f:
            for chunk in r.iter_content(chunk_size=1<<20):
                if chunk: f.write(chunk)
    log(f"✔ Salvo em: {path}")

def ensure_csi_index(vcf_path, force_reindex=False):
    csi = vcf_path + ".csi"
    tbi = vcf_path + ".tbi"
    if not force_reindex and (os.path.exists(csi) or os.path.exists(tbi)):
        log(f"✓ Índice presente para {os.path.basename(vcf_path)}")
        return
    log(f"🛠️  Gerando índice CSI para: {os.path.basename(vcf_path)}")
    pysam.tabix_index(vcf_path, preset="vcf", csi=True, force=True)
    log(f"✔ Índice criado: {csi}")

def pick_contig_alias(vcf_obj, target_chr: str):
    seqnames = list(vcf_obj.seqnames)
    candidates = [target_chr, f"chr{target_chr}", f"NC_{int(target_chr):06d}.11", f"NC_{int(target_chr):06d}.12"]
    for c in candidates:
        if c in seqnames:
            return c
    pat = re.compile(rf'(?:^|[^0-9]){re.escape(target_chr)}(?:[^0-9]|$)')
    for s in seqnames:
        if pat.search(s) and not any(x in s.lower() for x in ["random","decoy","unlocalized","unplaced"]):
            return s
    return target_chr

def normalize_clnsig(raw):
    if not raw:
        return None
    s = str(raw).replace("|", ",").replace(";", ",")
    vals = [v.strip() for v in s.split(",") if v.strip()]
    if not vals:
        return None
    out = []
    for v in vals:
        key = v.lower().replace("-", "_")
        out.append(CS_MAP.get(key, v.strip().title()))
    out = sorted(set(out))
    return json.dumps(out, ensure_ascii=False) if out else None

def infer_var_class(ref, alts):
    """Retorna 'SNP' para single-nucleotide; insertion/deletion/indel para demais."""
    if not ref or not alts:
        return "indel"
    ref = str(ref)
    alt_lens = [len(a) for a in alts if a and a != "."]
    if not alt_lens:
        return "indel"
    if len(ref) == 1 and all(l == 1 for l in alt_lens):
        return "SNP"
    gt = any(l > len(ref) for l in alt_lens)
    lt = any(l < len(ref) for l in alt_lens)
    if gt and not lt: return "insertion"
    if lt and not gt: return "deletion"
    return "indel"

def normalize_var_class(vc, ref, alts):
    """Normaliza CLNVC/valor livre -> classes curtas; single_nucleotide_variant/snv => SNP."""
    if vc:
        v = str(vc).strip().lower()
        if v in ("single_nucleotide_variant", "single nucleotide variant", "snv"):
            return "SNP"
        if v in ("insertion", "ins"):
            return "insertion"
        if v in ("deletion", "del"):
            return "deletion"
        if v in ("indel", "in-del", "in_del"):
            return "indel"
        # fallback pela estrutura:
        if ref and alts and len(ref) == 1 and all(len(a) == 1 for a in alts):
            return "SNP"
        return str(vc)
    return infer_var_class(ref, alts)

def parse_clinvar_mc(mc_str):
    """ClinVar INFO.MC -> JSON com consequências (ex.: 'missense_variant')."""
    if not mc_str:
        return None
    items = [x.strip() for x in str(mc_str).split(",") if x.strip()]
    labels = []
    for it in items:
        if "|" in it:
            _, name = it.split("|", 1)
            name = name.strip()
            if name:
                labels.append(name)
        else:
            labels.append(it.strip())
    if not labels:
        return None
    labels = sorted(set(labels))
    return json.dumps(labels, ensure_ascii=False)

def get_csq_index(vcf, field="Consequence"):
    raw = getattr(vcf, "raw_header", "")
    for line in raw.splitlines():
        if line.startswith("##INFO=<ID=CSQ"):
            if "Format:" in line:
                fmt = line.split("Format:")[-1].strip(" \">")
                fields = [f.strip() for f in fmt.split("|")]
                if field in fields:
                    return fields.index(field)
    return None

def parse_csq_all(csq_str, idx):
    if not csq_str or idx is None:
        return None
    consequences = []
    for item in str(csq_str).split(","):
        parts = item.split("|")
        if len(parts) > idx and parts[idx]:
            consequences.append(parts[idx])
    if not consequences:
        return None
    consequences = sorted(set(consequences))
    return json.dumps(consequences, ensure_ascii=False)

# ---------- Núcleo ----------

def run_pipeline(
    gene: str,
    out_sql_path: str,
    base_dir: str = "/data",
    skip_download: bool = False,
    force_reindex: bool = False,
    vcf_clinvar_path: str = None,
    vcf_ensembl_path: str = None,
):
    os.makedirs(base_dir, exist_ok=True)

    # 1) Coordenadas do gene (1 chamada REST)
    gene_url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}?content-type=application/json"
    log(f"[Ensembl] GET {gene_url}")
    r = requests.get(gene_url, headers={"Content-Type":"application/json"})
    r.raise_for_status()
    gd = r.json()
    chr_   = str(gd["seq_region_name"])
    start  = int(gd["start"])
    end    = int(gd["end"])
    strand = int(gd.get("strand", 0))
    log(f"Gene {gene}: chr{chr_}:{start}-{end} | strand={strand}")

    # 2) Caminhos dos VCFs
    if not vcf_clinvar_path:
        vcf_clinvar_path = os.path.join(base_dir, "clinvar.vcf.gz")
        clinvar_url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
        download_if_needed(clinvar_url, vcf_clinvar_path, skip_download=skip_download)
    else:
        log(f"Usando ClinVar VCF: {vcf_clinvar_path}")

    if not vcf_ensembl_path:
        vcf_ensembl_path = os.path.join(base_dir, f"homo_sapiens-chr{chr_}.vcf.gz")
        ensembl_url = f"https://ftp.ensembl.org/pub/release-113/variation/vcf/homo_sapiens/homo_sapiens-chr{chr_}.vcf.gz"
        download_if_needed(ensembl_url, vcf_ensembl_path, skip_download=skip_download)
    else:
        log(f"Usando Ensembl VCF: {vcf_ensembl_path}")

    # 3) Indexação CSI (se necessário)
    ensure_csi_index(vcf_clinvar_path, force_reindex=force_reindex)
    ensure_csi_index(vcf_ensembl_path, force_reindex=force_reindex)

    # 4) ClinVar: extrair região
    vcf_cln = VCF(vcf_clinvar_path)
    contig_cln = pick_contig_alias(vcf_cln, chr_)
    region_cln = f"{contig_cln}:{start}-{end}"
    log(f"[ClinVar] região: {region_cln}")

    clin_rows = []
    for rec in vcf_cln(region_cln):
        info = rec.INFO or {}
        rsid = rec.ID or ""
        var_class = normalize_var_class(info.get("CLNVC"), rec.REF, rec.ALT or [])
        clnsig_json = normalize_clnsig(info.get("CLNSIG"))
        cons_json   = parse_clinvar_mc(info.get("MC"))
        clin_rows.append({
            "key": (contig_cln, int(rec.POS), rec.REF, tuple(rec.ALT or [])),
            "variant_id": rsid,                         # pode estar vazio
            "gene_display_name": gene,
            "consequence_type": cons_json or "",
            "start_pos": int(rec.POS),
            "end_pos": int(rec.end),
            "strand": strand if strand in (-1, 1) else None,
            "seq_region_name": contig_cln,
            "assembly_name": "GRCh38",
            "source": "ClinVar",
            "var_class": var_class,
            "clinical_significance": clnsig_json or ""
        })
    df_cln = pd.DataFrame(clin_rows)

    # 5) Ensembl: fallback de consequence
    vcf_ens = VCF(vcf_ensembl_path)
    contig_ens = pick_contig_alias(vcf_ens, chr_)
    region_ens = f"{contig_ens}:{start}-{end}"
    log(f"[Ensembl] região: {region_ens}")

    csq_idx = get_csq_index(vcf_ens, "Consequence")
    ens_map = {}
    for rec in vcf_ens(region_ens):
        cons_json = parse_csq_all(rec.INFO.get("CSQ"), csq_idx) if "CSQ" in rec.INFO else None
        ens_map[(contig_ens, int(rec.POS), rec.REF, tuple(rec.ALT or []))] = cons_json or ""

    # 6) Integrar e filtrar por clinical_significance
    final_rows = []
    for r in clin_rows:
        if not r["clinical_significance"]:
            continue
        if not r["consequence_type"]:
            k2 = (contig_ens, r["start_pos"], r["key"][2], r["key"][3])
            r["consequence_type"] = ens_map.get(k2, "")
        final_rows.append(r)

    final = pd.DataFrame(final_rows).fillna("")
    log(f"Final (com clinical_significance): {len(final)}")

    # 7) Gerar SQL único
    lines = []
    lines.append("-- Variants + Variant Alleles (único arquivo)")
    lines.append("START TRANSACTION;")
    lines.append("SET NAMES utf8mb4;")
    lines.append("SET time_zone = '+00:00';")
    lines.append("SET FOREIGN_KEY_CHECKS=0;")

    def variant_insert_sql(r):
        variant_id = r["variant_id"] or make_synth_id(r["key"])
        cols = {
            "variant_id":        f"'{sql_escape(variant_id)}'",
            "gene_display_name": f"'{sql_escape(r['gene_display_name'])}'",
            "consequence_type":  f"'{sql_escape(r['consequence_type'])}'" if r["consequence_type"] else "NULL",
            "start_pos":         str(int(r["start_pos"])) if str(r["start_pos"]).isdigit() else "NULL",
            "end_pos":           str(int(r["end_pos"]))   if str(r["end_pos"]).isdigit()   else "NULL",
            "strand":            str(int(r["strand"]))    if str(r["strand"]) in ("-1","1") else "NULL",
            "seq_region_name":   f"'{sql_escape(r['seq_region_name'])}'" if r["seq_region_name"] else "NULL",
            "assembly_name":     f"'{sql_escape(r['assembly_name'])}'"   if r["assembly_name"]   else "NULL",
            "source":            f"'{sql_escape(r['source'])}'"          if r["source"]          else "NULL",
            "clinical_significance": f"CAST('{sql_escape(r['clinical_significance'])}' AS JSON)" if r["clinical_significance"] else "NULL",
            "var_class":         f"'{sql_escape(r['var_class'])}'" if r["var_class"] else "NULL",
        }
        col_names = ", ".join(cols.keys())
        col_vals  = ", ".join(cols.values())
        upd = (
            "consequence_type=VALUES(consequence_type), "
            "start_pos=VALUES(start_pos), end_pos=VALUES(end_pos), strand=VALUES(strand), "
            "seq_region_name=VALUES(seq_region_name), assembly_name=VALUES(assembly_name), "
            "source=VALUES(source), clinical_significance=VALUES(clinical_significance), "
            "var_class=VALUES(var_class)"
        )
        return f"INSERT INTO variant ({col_names}) VALUES ({col_vals}) ON DUPLICATE KEY UPDATE {upd};"

    def variant_allele_inserts_sql(r):
        inserts = []
        vid = r["variant_id"] or make_synth_id(r["key"])
        gene = r["gene_display_name"]
        if not isinstance(r["key"], tuple) or len(r["key"]) != 4:
            return inserts
        contig, pos, ref, alts = r["key"]
        if ref and ref != ".":
            inserts.append(
                "INSERT INTO variant_allele (variant_id, gene_display_name, allele, is_reference) "
                f"VALUES ('{sql_escape(vid)}','{sql_escape(gene)}','{sql_escape(ref)}',1) "
                "ON DUPLICATE KEY UPDATE is_reference=VALUES(is_reference);"
            )
        for a in (alts or []):
            if a and a != ".":
                inserts.append(
                    "INSERT INTO variant_allele (variant_id, gene_display_name, allele, is_reference) "
                    f"VALUES ('{sql_escape(vid)}','{sql_escape(gene)}','{sql_escape(a)}',0) "
                    "ON DUPLICATE KEY UPDATE is_reference=VALUES(is_reference);"
                )
        return inserts

    for _, r in final.iterrows():
        lines.append(variant_insert_sql(r))

    seen_alleles = set()
    for _, r in final.iterrows():
        for stmt in variant_allele_inserts_sql(r):
            key = tuple(re.findall(r"VALUES \('([^']*)','([^']*)','([^']*)',([01])\)", stmt)[0])
            if key not in seen_alleles:
                seen_alleles.add(key)
                lines.append(stmt)

    lines.append("SET FOREIGN_KEY_CHECKS=1;")
    lines.append("COMMIT;")

    with open(out_sql_path, "w") as f:
        f.write("\n".join(lines))

    log(f"✔ SQL único salvo em: {out_sql_path}")
    log(f"Linhas no SQL: {len(lines)}")

# ---------- CLI ----------

def parse_args():
    p = argparse.ArgumentParser(description="ClinVar+Ensembl -> SQL (variant + variant_allele)")
    p.add_argument("--gene", "-g", required=True, help="Símbolo do gene (ex.: BRCA2)")
    p.add_argument("--out-sql", default="/data/variants_full.sql", help="Caminho de saída do .sql")
    p.add_argument("--base-dir", default="/data", help="Pasta base para VCFs/índices")
    p.add_argument("--skip-download", action="store_true", help="Não baixar arquivos se não existirem (erro se ausentes)")
    p.add_argument("--force-reindex", action="store_true", help="Sempre recriar .csi")
    p.add_argument("--vcf-clinvar", help="Path para clinvar.vcf.gz (se já tiver local)")
    p.add_argument("--vcf-ensembl", help="Path para homo_sapiens-chrN.vcf.gz (se já tiver local)")
    return p.parse_args()

if __name__ == "__main__":
    args = parse_args()
    try:
        run_pipeline(
            gene=args.gene,
            out_sql_path=args.out_sql,
            base_dir=args.base_dir,
            skip_download=args.skip_download,
            force_reindex=args.force_reindex,
            vcf_clinvar_path=args.vcf_clinvar,
            vcf_ensembl_path=args.vcf_ensembl,
        )
    except Exception as e:
        log(f"[ERRO] {e}")
        sys.exit(1)