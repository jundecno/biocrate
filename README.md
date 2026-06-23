# biocrate

`biocrate` is a small Python toolkit for frequently used bioinformatics and
molecular data operations, including FASTA/JSON conversion, structure format
conversion, pocket residue extraction, simple file IO helpers, and common model
evaluation metrics.

The package uses a `src/` layout and keeps the top-level import lightweight:

```python
import biocrate
from biocrate.metrics import evaluate_binary_classification
from biocrate.general import fasta_load, json_dump
```

## Install

For local development:

```bash
pip install -e ".[dev]"
```

For chemistry-related helpers:

```bash
pip install -e ".[chem]"
```

For OpenBabel-based format conversion:

```bash
pip install -e ".[openbabel]"
```

## Modules

- `biocrate.general`: file IO, FASTA helpers, LMDB helpers, structure load/save
- `biocrate.convertors`: FASTA, PDB/mmCIF, SDF/MOL2 and sequence format conversion
- `biocrate.operators`: structure and path manipulation helpers
- `biocrate.calculators`: pocket residue extraction and molecule conformation generation
- `biocrate.fetch`: lightweight download helpers for public biological resources
- `biocrate.metrics`: common classification and regression metrics
- `biocrate.constants`: shared constants and mirrors

## Test

```bash
pytest
```
