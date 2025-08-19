FROM python:3.11-slim

# deps nativas para cyvcf2/pysam
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev ca-certificates \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /app
ENV BASE_DIR=/data

# dependências
COPY requirements.txt /app/
RUN pip install --no-cache-dir -r requirements.txt

# código
COPY pipeline.py /app/pipeline.py
COPY app.py /app/app.py

# diretório de dados (montável) e porta
VOLUME ["/data"]
EXPOSE 8000

# rode FastAPI (Uvicorn)
CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "8000", "--timeout-keep-alive", "120"]