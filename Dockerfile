FROM python:3.11-slim
WORKDIR /work
COPY pyproject.toml README.md ./
COPY src ./src
RUN python -m pip install --no-cache-dir .
ENTRYPOINT ["pbmc3k"]
