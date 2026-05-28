# Batteries-included image for ISCN-parser.
# Build:  docker build -t iscn-parser .
# Run:    docker run --rm iscn-parser --string "arr[GRCh37] 1p36.33p36.31(849466_6015936)x1"
FROM python:3.12-slim

LABEL org.opencontainers.image.title="ISCN-parser" \
      org.opencontainers.image.description="Parse ISCN cytogenomic nomenclature and convert to PLINK .cnv or BED." \
      org.opencontainers.image.source="https://github.com/jlanej/ISCN-parser" \
      org.opencontainers.image.licenses="MIT"

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /app

# Install the package. Copy metadata first to leverage Docker layer caching.
COPY pyproject.toml README.md LICENSE ./
COPY src ./src
RUN pip install --no-cache-dir .

# Bundle example data so the image is usable out of the box.
COPY examples ./examples

# Run as a non-root user.
RUN useradd --create-home --uid 1000 iscn
USER iscn

ENTRYPOINT ["iscn-parser"]
CMD ["--help"]
