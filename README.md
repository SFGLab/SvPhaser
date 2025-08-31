# SvPhaser

> **Haplotype‑aware structural‑variant phasing for long‑read data**

[![PyPI version](https://img.shields.io/pypi/v/svphaser.svg?logo=pypi)](https://pypi.org/project/svphaser)
[![Tests](https://img.shields.io/github/actions/workflow/status/your-org/SvPhaser/ci.yml?label=ci)](https://github.com/your-org/SvPhaser/actions)
[![License](https://img.shields.io/github/license/your-org/SvPhaser.svg)](LICENSE)

---

`SvPhaser` assigns **haplotype genotypes** to *pre‑called structural variants (SVs)* using *HP‑tagged* long‑read alignments (ONT, PacBio). Think of it as *WhatsHap‑style phasing* for DEL/INS/INV/…: we don’t discover SVs; we **phase** them — outputting `1|0`, `0|1`, `1|1`, or `./.` plus a **Genotype Quality (GQ)**.

## Highlights

* **Fast, per‑chromosome multiprocessing** with linear scale‑out.
* **Deterministic decision tree** with clear thresholds (`major_delta`, `equal_delta`).
* **Overflow‑safe GQ** (exact binomial for shallow depth; normal approx for deep coverage; **capped at 99**).
* **Standards‑friendly VCF**: preserves your original header/INFO; writes `GT:GQ` in **FORMAT** and optional `GQBIN` in **INFO**; emits clean, tab‑delimited records.
* **Transparent QC**: CSV summary, optional GQ bin labels, and a `_dropped_svs.csv` audit file.
* **Simple CLI & Python API** with sane defaults.

---

## Installation

```bash
# Requires Python ≥3.9
pip install svphaser              # from PyPI
# or a specific tag from source
pip install "git+https://github.com/your-org/SvPhaser.git@v2.0.1"
```

Runtime deps (`cyvcf2`, `pysam`, `pandas`, `typer`) are installed automatically.

---

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

Outputs written to **`results/`**:

```
sample_unphased_phased.vcf        # original VCF + GT:GQ (FORMAT) and optional GQBIN (INFO)
sample_unphased_phased.csv        # tidy table: chrom,pos,id,svtype,n1,n2,gt,gq,gq_label
sample_unphased_dropped_svs.csv   # variants removed by global depth filter (for audit)
```

### Genotype model (at a glance)

For each SV, count HP‑tagged reads across its span → `n1` (HP1), `n2` (HP2), `N=n1+n2`.

* **Global support filter**: drop only if `(n1 < min_support) AND (n2 < min_support)`.
* **Decision tree** (ratios):

  * `n1/N ≥ major_delta` → `1|0`
  * `n2/N ≥ major_delta` → `0|1`
  * `|n1−n2|/N ≤ equal_delta` → `1|1`
  * else → `./.`
* **Genotype quality**: exact binomial tail for `N ≤ 200`; continuity‑corrected normal approx for `N > 200`; Phred‑scaled, **capped at 99**. Optional `GQBIN` label via `--gq-bins` (e.g., `"30:High,10:Moderate"`).

Default knobs: `--min-support 10`, `--major-delta 0.70`, `--equal-delta 0.25`.

---

## Python API

Two entry points — a convenience wrapper and the full engine:

```python
from pathlib import Path
from svphaser import phase  # convenience wrapper

out_vcf, out_csv = phase(
    sv_vcf="sample.vcf.gz",
    bam="sample.bam",
    out_dir="results",
    min_support=10,
    major_delta=0.70,
    equal_delta=0.25,
    gq_bins="30:High,10:Moderate",
    threads=8,
)
```

Or call the high‑level engine explicitly:

```python
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

---

## Repo layout

```
SvPhaser/
├─ src/svphaser/        # package
│  ├─ cli.py            # Typer entry‑point (svphaser phase …)
│  ├─ logging.py        # minimal logging setup
│  └─ phasing/
│     ├─ algorithms.py  # decision tree + overflow‑safe GQ
│     ├─ io.py          # driver & VCF/CSV writer
│     ├─ _workers.py    # per‑chromosome workers
│     └─ types.py       # dataclasses & aliases
├─ tests/               # pytest + smoke tests
├─ docs/                # extra documentation
└─ docs/result_images/  # generated plots & diagrams
```

---

## Development

```bash
git clone https://github.com/your-org/SvPhaser.git && cd SvPhaser
python -m venv .venv && source .venv/bin/activate
pip install -e .[dev]

# one‑time
pre-commit install

# checks
pre-commit run --all-files
pytest -q
mypy src/svphaser
```

We use **hatch‑vcs** for versioning — releases are driven by Git tags like `v2.0.1`.

---

## Citing SvPhaser

If SvPhaser contributed to your research, please cite:

```bibtex
@software{svphaser2025,
  author       = {Pranjul Mishra and Sachin Gadakh},
  title        = {SvPhaser: haplotype‑aware structural‑variant phasing},
  version      = {2.0.1},
  date         = {2025-08-31},
  url          = {https://github.com/your-org/SvPhaser}
}
```

---

## License

`SvPhaser` is released under the MIT License – see [LICENSE](LICENSE).

## 📬 Contact

Developed by **Team5** (*BioAI Hackathon*) – Sachin Gadakh & Pranjul Mishra.

Lead contacts:

* [pranjul.mishra@proton.me](mailto:pranjul.mishra@proton.me)
* [s.gadakh@cent.uw.edu.pl](mailto:s.gadakh@cent.uw.edu.pl)

Feedback, feature requests, and bug reports are welcome — please open a GitHub issue or reach out by e‑mail.

---

*Happy phasing!*
