# app.py
import os
import time
from pathlib import Path
from fastapi import FastAPI, Query, HTTPException
from fastapi.responses import FileResponse, JSONResponse
from fastapi.middleware.cors import CORSMiddleware

# importe seu pipeline
from pipeline import run_pipeline

app = FastAPI(title="Clinical Significant Variants from Gene")

# CORS (ajuste os domínios conforme precisar)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"], allow_credentials=True,
    allow_methods=["*"], allow_headers=["*"],
)

DATA_DIR = Path(os.environ.get("BASE_DIR", "/data"))

@app.get("/health")
def health():
    return {"status": "ok"}

@app.get("/run")
def run(
    gene: str = Query(..., min_length=2, description="Símbolo do gene, ex.: BRCA2"),
    skip_download: bool = Query(False, description="Pular downloads se VCFs/índices já existirem"),
):
    """
    Executa o pipeline e retorna o arquivo .sql como download.
    """
    try:
        DATA_DIR.mkdir(parents=True, exist_ok=True)
        ts = time.strftime("%Y%m%d-%H%M%S")
        out_sql = DATA_DIR / f"variants_full_{gene}_{ts}.sql"

        run_pipeline(
            gene=gene,
            out_sql_path=str(out_sql),
            base_dir=str(DATA_DIR),
            skip_download=skip_download,
            force_reindex=False,
            vcf_clinvar_path=None,
            vcf_ensembl_path=None,
        )

        if not out_sql.exists():
            raise HTTPException(status_code=500, detail="Falha ao gerar SQL")

        # retorna para download
        return FileResponse(
            path=str(out_sql),
            filename=out_sql.name,
            media_type="application/sql",
        )
    except Exception as e:
        return JSONResponse(status_code=500, content={"error": str(e)})