# SvPhaser

> **Haplotype-aware structural-variant (SV) genotyper for long-read data**

[![PyPI version](https://img.shields.io/pypi/v/svphaser.svg?logo=pypi)](https://pypi.org/project/svphaser/)
[![Python](https://img.shields.io/pypi/pyversions/svphaser.svg)](https://pypi.org/project/svphaser/)
[![License](https://img.shields.io/github/license/SFGLab/SvPhaser.svg)](LICENSE)

---

**SvPhaser** phases **pre-called structural variants (SVs)** using **HP-tagged** long-read alignments (PacBio HiFi, ONT Q20+, …).

Think of it as *WhatsHap* for insertions/deletions/duplications:
- **we do not discover SVs**
- **we assign haplotype genotypes** (`0|1`, `1|0`, `1|1`, or `./.`)
- and compute a **Genotype Quality (GQ)** score

All in a single, embarrassingly-parallel pass over the genome.

## Highlights

- **Fast per-chromosome multiprocessing** (scale-out on multi-core CPUs).
- **Deterministic Δ-based decision logic** (no MCMC / HMM).
- **CLI + Python API**.
- **Non-destructive VCF augmentation**: injects phasing fields while preserving the original header and records.
- **Configurable confidence bins** + optional plots.

## Installation

### From PyPI (recommended)

```bash
# Requires Python >= 3.9
pip install svphaser
````

Optional extras (if you use them):

```bash
pip install "svphaser[plots]"
pip install "svphaser[bench]"
pip install "svphaser[dev]"
```

### From source

```bash
git clone https://github.com/SFGLab/SvPhaser.git
cd SvPhaser
pip install -e .
```

## Inputs & requirements

SvPhaser expects:

1. **Unphased SV VCF** (`.vcf` / `.vcf.gz`)

   * SVs should already be called by your preferred SV caller.

2. **HP-tagged BAM** (long-read alignments)

   * Reads must contain haplotype tags (e.g., `HP`) produced by an upstream phasing pipeline.

If your BAM is not HP-tagged, SvPhaser cannot assign haplotypes.

## Quick start (CLI)

```bash
svphaser phase \
    sample_unphased.vcf.gz \
    sample.sorted_phased.bam \
    --out-dir results/ \
    --min-support 10 \
    --major-delta 0.70 \
    --equal-delta 0.25 \
    --gq-bins "30:High,10:Moderate" \
    --threads 32
```

### Outputs

Inside `results/`:

* `*_phased.vcf` — your original VCF with additional INFO fields:

  * `HP_GT` — phased genotype
  * `HP_GQ` — genotype quality score
  * `HP_GQBIN` — confidence bin label (based on your `--gq-bins`)
* `*_phased.csv` — tidy table for plotting / downstream analysis

For algorithmic details, see: **`docs/methodology.md`**.

## Python API

```python
from pathlib import Path
from svphaser.phasing.io import phase_vcf

phase_vcf(
    Path("sample.vcf.gz"),
    Path("sample.bam"),
    out_dir=Path("results"),
    min_support=10,
    major_delta=0.70,
    equal_delta=0.25,
    gq_bins="30:High,10:Moderate",
    threads=8,
)
```

The phased table can also be loaded from the generated CSV for custom analytics.

## Repository structure (high level)

```
SvPhaser/
├─ src/svphaser/         # importable package
├─ tests/                # test suite + small fixtures (if present)
├─ docs/                 # methodology + notes
├─ notebooks/            # experiments / analysis (if present)
├─ figures/              # plots & diagrams (if present)
├─ pyproject.toml
└─ CHANGELOG.md
```

## Development

```bash
git clone https://github.com/SFGLab/SvPhaser.git
cd SvPhaser

python -m venv .venv
source .venv/bin/activate

pip install -e ".[dev]"
pytest -q
mypy src/svphaser
```

See `CONTRIBUTING.md` for contribution guidelines.

## Citing SvPhaser

If SvPhaser contributed to your research, please cite:

```bibtex
@software{svphaser2025,
  author  = {Pranjul Mishra and Sachin Gadakh},
  title   = {SvPhaser: Haplotype-aware structural-variant genotyping from HP-tagged long-read BAMs},
  version = {2.0.6},
  year    = {2025},
  month   = nov,
  url     = {https://github.com/SFGLab/SvPhaser},
  note    = {PyPI: https://pypi.org/project/svphaser/}
}
```

(If you need maximum rigor for a paper, cite a specific git commit hash too.)

## License

SvPhaser is released under the **MIT License** — see [LICENSE](LICENSE).

## Contact

Developed by **Team 5 (BioAI Hackathon)**.

* Pranjul Mishra — [pranjul.mishra@proton.me](mailto:pranjul.mishra@proton.me)
* Sachin Gadakh — [s.gadakh@cent.uw.edu.pl](mailto:s.gadakh@cent.uw.edu.pl)

Issues and feature requests: please open a GitHub issue.

```

### Two hard notes (don’t ignore)
- If you **don’t actually have CI set up**, don’t show a CI badge. A fake badge is worse than no badge.
- If your repo layout doesn’t include `notebooks/figures/tests fixtures`, either adjust that tree block or remove it to avoid “template smell.”

If you want, paste your **current `.github/workflows` filenames** (or tell me if you have none) and I’ll add the *correct* CI badge line too—without guessing.
::contentReference[oaicite:1]{index=1}
```

[1]: https://pypi.org/project/svphaser/ "svphaser · PyPI"
