# ISCN-parser

Parse [**ISCN**](https://iscn.karger.com/) (the *International System for Human
Cytogenomic Nomenclature*) results and convert them into
[PLINK `.cnv`](https://www.cog-genomics.org/plink/1.9/formats#cnv) segment files
or [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) intervals.

ISCN is the standard notation cytogeneticists use to describe karyotypes and
microarray results (see the
[Wikipedia overview](https://en.wikipedia.org/wiki/International_System_for_Human_Cytogenomic_Nomenclature)).
This tool turns those human-readable strings into the tabular,
coordinate-based formats used by downstream genomics pipelines.

[![CI](https://github.com/jlanej/ISCN-parser/actions/workflows/ci.yml/badge.svg)](https://github.com/jlanej/ISCN-parser/actions/workflows/ci.yml)

---

## Features

- **Microarray (`arr`) nomenclature** parsed losslessly to coordinates, e.g.
  `arr[GRCh37] 1p36.33p36.31(849466_6015936)x1`. Handles:
  - genome build detection (`[GRCh37]`, `[hg19]`, `[GRCh38]`, `[hg38]`, `[NCBI36]`/`[hg18]`);
  - coordinate separators `_` (current) and `-` (legacy), with comma thousands separators;
  - copy numbers `x0`–`x9` (homozygous/heterozygous loss, gains);
  - mosaic ranges (`x1~2`) and mosaic fractions (`x1[0.45]`);
  - copy-neutral loss of heterozygosity (`... x2 hmz`);
  - multiple comma-separated aberrations in one string.
- **Conventional karyotypes**, e.g. `47,XX,+21`, `45,X`, `47,XXY`,
  `46,XY,del(5)(p15.2p15.33)`:
  - whole-chromosome aneuploidies (`+21`, `-7`) resolved to full chromosome spans;
  - sex-chromosome complement changes (Turner, Klinefelter, XYY, XXX, …);
  - `del`/`dup` band events resolved with the bundled UCSC cytoBand files
    (or a custom file via `--cytoband`; falls back to imprecise arm-level
    coordinates with a warning when no file is provided);
  - balanced rearrangements (`t`, `inv`, `ins`) recognised as copy-number neutral.
- **Robust, non-fatal error handling**: malformed records yield diagnostics, not
  crashes, so a single bad line never aborts a whole file.
- **Outputs**: PLINK `.cnv` (+ companion `.fam` and `.cnv.map`) and colour-coded BED9.
- **Batteries included**: bundled reference data (chromosome lengths, centromere
  positions, and full UCSC cytoBand tables for GRCh37 and GRCh38), example
  datasets, a CLI, a Python API, a Docker image, and CI covering all the
  variations above.

## Installation

```bash
pip install .            # from a clone of this repository
# or, for development:
pip install -e ".[dev]"
```

Requires Python 3.9+. There are no runtime dependencies.

### Docker

A batteries-included image is built and published to the GitHub Container
Registry by CI:

```bash
docker pull ghcr.io/jlanej/iscn-parser:latest
docker run --rm ghcr.io/jlanej/iscn-parser:latest \
  --string "arr[GRCh37] 1p36.33p36.31(849466_6015936)x1"
```

To build it locally:

```bash
docker build -t iscn-parser .
docker run --rm iscn-parser --help
```

## Usage

### Command line

Convert a single string to PLINK CNV (written to stdout):

```bash
iscn-parser --string "arr[GRCh37] 1p36.33p36.31(849466_6015936)x1"
```

Convert a file of records to PLINK output files (`<prefix>.cnv`, `.fam`, `.cnv.map`):

```bash
iscn-parser --input examples/array_samples.tsv --format cnv --output results/array
```

Convert to BED:

```bash
iscn-parser --input examples/array_samples.tsv --format bed --output results/array.bed
```

Convert conventional karyotypes (a genome build is required to resolve coordinates):

```bash
iscn-parser --input examples/karyotype_samples.tsv --build GRCh37 --output results/kar
```

Resolve band-level `del`/`dup` precisely using the bundled cytoBand files:

```bash
# From a source checkout or editable install:
iscn-parser --string "46,XY,del(5)(p15.2p15.33)" --build GRCh37 \
  --cytoband src/iscn_parser/data/cytoBand_GRCh37.txt.gz

# From any install, resolve the bundled file path via Python:
CYTOBAND=$(python -c "from importlib.resources import files; \
  print(files('iscn_parser').joinpath('data/cytoBand_GRCh37.txt.gz'))")
iscn-parser --string "46,XY,del(5)(p15.2p15.33)" --build GRCh37 \
  --cytoband "$CYTOBAND"
```

Or supply a custom/updated cytoBand file:

```bash
iscn-parser --string "46,XY,del(5)(p15.2p15.33)" --build GRCh37 \
  --cytoband cytoBand.txt.gz
```

Useful flags: `--strict` (exit non-zero on parse errors), `--quiet` (suppress
diagnostics), `--sample-column`/`--iscn-column`/`--delimiter`/`--has-header`
(control input layout), `--no-map` (skip the `.cnv.map` file). Run
`iscn-parser --help` for the full list.

#### Input file format

One record per line. A line may be either a bare ISCN string or a delimited
record carrying a sample id (default: tab-separated, sample id in column 0,
ISCN string in column 1). Blank lines and `#` comments are ignored.

```
ARR001	arr[GRCh37] 1p36.33p36.31(849466_6015936)x1
ARR002	arr[GRCh37] 7q11.23(72726578_74139390)x3
```

### Python API

```python
from iscn_parser import parse_iscn, to_plink_cnv, to_bed, GenomeBuild

result = parse_iscn("arr[GRCh37] 1p36.33p36.31(849466_6015936)x1")
for v in result.variants:
    print(v.chrom, v.start, v.end, v.copy_number, v.type_label)

print(to_plink_cnv(result.variants))
print(to_bed(result.variants))

# Files of records, with a default build for karyotypes:
from iscn_parser import parse_iscn_file, write_plink
result = parse_iscn_file("examples/karyotype_samples.tsv", build=GenomeBuild.GRCh37)
write_plink(result.variants, "results/kar")
```

`result.messages` contains `INFO`/`WARNING`/`ERROR` diagnostics; `result.ok` is
`True` when no errors occurred.

## Output formats

### PLINK `.cnv`

A whitespace-delimited segment file with the columns
`FID IID CHR BP1 BP2 TYPE SCORE SITES`. Following the PLINK convention, **`TYPE`
holds the integer copy number** (`0` homozygous deletion, `1` heterozygous
deletion, `2` diploid, `3`/`4`… gains). Coordinates are 1-based inclusive and
the sex chromosomes are recoded numerically (`X→23`, `Y→24`, `XY→25`, `MT→26`).
Companion `.fam` (sample list) and `.cnv.map` (breakpoint markers) files are
written alongside it. See the
[PLINK CNV documentation](https://www.cog-genomics.org/plink/1.9/formats#cnv).

### BED

A 9-column BED file. ISCN's 1-based inclusive coordinates are converted to BED's
0-based half-open convention. The `name` column encodes the event type and copy
number (e.g. `DEL|CN1|1p36.33p36.31`), and `itemRgb` colours losses red, gains
blue, and copy-neutral LOH green for genome-browser display.

## Methods & how parsing works

1. **Build detection.** An embedded build token (`[GRCh37]`, `[hg19]`, …) takes
   precedence over the `--build` default. Aliases are normalised to canonical
   builds (GRCh37/GRCh38/NCBI36).
2. **Tokenisation.** A string is split on commas *only at bracket depth zero*, so
   coordinate groups containing comma thousands separators
   (`(849,466_6,015,936)`) stay intact.
3. **Dispatch.** Strings beginning with `arr` are parsed as microarray results;
   everything else is treated as a conventional karyotype (modal number, sex
   complement, then aberrations).
4. **Array aberrations** are matched by a tolerant regular expression capturing
   chromosome, optional bands, base-pair start/end, copy number, optional mosaic
   range/fraction, and trailing modifiers (`hmz`, `loh`). Copy-neutral `x2`
   regions are skipped unless flagged as LOH.
5. **Karyotype aberrations**: `+N`/`-N` become whole-chromosome gains/losses
   using bundled chromosome lengths; sex-complement changes are computed against
   the nearest normal complement; `del`/`dup` bands are resolved using the
   bundled cytoBand table for the detected build (or a custom file supplied via
   `--cytoband`), falling back to an arm-level approximation with a warning
   when no matching band is found. Balanced rearrangements are reported as
   copy-number neutral.
6. **Coordinate model.** Internally all events are 1-based inclusive; BED output
   converts to 0-based half-open at the boundary.

### Reference data

Chromosome lengths, approximate centromere boundaries, and full **UCSC cytoBand
tables** for **GRCh37** and **GRCh38** are bundled in `src/iscn_parser/data/`:

| File | Purpose |
|------|---------|
| `chromosomes_GRCh37.tsv` | Chromosome lengths + centromere positions (GRCh37) |
| `chromosomes_GRCh38.tsv` | Chromosome lengths + centromere positions (GRCh38) |
| `cytoBand_GRCh37.txt.gz` | Full UCSC cytoBand table for GRCh37/hg19 (862 bands) |
| `cytoBand_GRCh38.txt.gz` | Full UCSC cytoBand table for GRCh38/hg38 (1549 bands) |

Whole-chromosome and arm-level events are resolved from the chromosome tables
alone. The bundled cytoBand files enable **exact band-level coordinate
resolution** for `del`/`dup` karyotype events without any external downloads.
You can also supply a custom or updated cytoBand file via `--cytoband` if
needed; the script used to download the originals from UCSC is available at
`scripts/download_cytobands.py`.

## Limitations

- Band-level conventional events fall back to whole-arm coordinates (with a
  warning) only when neither the bundled cytoBand files nor a custom file can
  resolve the requested bands.
- Balanced translocations/inversions carry no copy-number change and so produce
  no CNV records by design.
- Complex constructs (derivative chromosomes, marker chromosomes, ring
  chromosomes, isochromosomes) are reported as unsupported rather than guessed.

## Development

```bash
pip install -e ".[dev]"
ruff check src tests
pytest
```

The test suite (`tests/`) exercises array, karyotype, malformed, converter,
CLI, reference-data, and **full cytoband** behaviour (loading, band/range
lookups, integration with the bundled GRCh37/GRCh38 cytoBand files, and
CLI `--cytoband` flag), and ships extensive example data under
`tests/data/` and `examples/`.

## License

[MIT](LICENSE).
