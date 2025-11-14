# SvPhaser

> **Haplotype‑aware structural‑variant genotyper for long‑read data**

[![PyPI version](https://img.shields.io/pypi/v/svphaser.svg?logo=pypi)](https://pypi.org/project/svphaser)
[![Tests](https://img.shields.io/github/actions/workflow/status/your‑org/SvPhaser/ci.yml?label=ci)](https://github.com/your‑org/SvPhaser/actions)
[![License](https://img.shields.io/github/license/your‑org/SvPhaser.svg)](LICENSE)

---

`SvPhaser` phases **pre‑called structural variants (SVs)** using *HP‑tagged* long‑read alignments (PacBio HiFi, ONT Q20+, …).  Think of it as *WhatsHap* for insertions/deletions/duplications: we do **not** discover SVs; we assign each variant a haplotype genotype (`0|1`, `1|0`, `1|1`, or `./.`) together with a **Genotype Quality (GQ)** score – all in a single, embarrassingly‑parallel pass over the genome.

## Key highlights

* **Fast, per‑chromosome multiprocessing** – linear scale‑out on 32‑core workstations.
* **Deterministic Δ‑based decision tree** – no MCMC or hidden state machines.
* **Friendly CLI** (`svphaser phase …`) and importable Python API.
* **Seamless VCF injection** – adds `HP_GT`, `HP_GQ`, `HP_GQBIN` INFO tags while copying the original header verbatim.
* **Configurable confidence bins** and publication‑ready plots (see `result_images/`).

---

## Installation

```bash
# Requires Python ≥3.9
pip install svphaser            # PyPI (coming soon)
# or
pip install git+https://github.com/your‑org/SvPhaser.git@v0.2.0
```

`cyvcf2`, `pysam`, `typer[all]`, and `pandas` are pulled in automatically.

## Quick‑start

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

Outputs (written inside **`results/`**)

```
sample_unphased_phased.vcf   # original VCF + HP_* INFO fields
sample_unphased_phased.csv   # tidy table for plotting / downstream R
```

See [`docs/methodology.md`](docs/Methodology.md) and the flow‑chart below for algorithmic details.

![SvPhaser methodology](docs/result_images/methodology_diagram.png)

## Folder layout

```
SvPhaser/
├─ src/svphaser/        # importable package
│  ├─ cli.py            # Typer entry‑point
│  ├─ logging.py        # unified log setup
│  └─ phasing/
│     ├─ algorithms.py  # core maths
│     ├─ io.py          # driver & I/O
│     ├─ _workers.py    # per‑chrom processes
│     └─ types.py       # thin dataclasses
├─ tests/               # pytest suite + mini data
├─ docs/                # extra documentation
├─ result_images/       # generated plots & diagrams
└─ CHANGELOG.md
```

## Python usage

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

The resulting `DataFrame` can be loaded from the CSV for custom analytics.




## Development & contributing

1. Clone and create a virtual env:

   ```bash
   git clone https://github.com/your‑org/SvPhaser.git && cd SvPhaser
   python -m venv .venv && source .venv/bin/activate
   pip install -e .[dev]
   ```
2. Run the test‑suite & type checks:

   ```bash
   pytest -q
   mypy src/svphaser
   black --check src tests
   ```
3. Send a PR targeting the **`dev`** branch; one topic per PR.

Please read `CONTRIBUTING.md` (to come) for style‑guides and the DCO sign‑off.

## Citing SvPhaser

If SvPhaser contributed to your research, please cite:

```bibtex
@software{svphaser2024,
  author       = {Pranjul Mishra, Sachin Ghadak, CeNT Lab},
  title        = {SvPhaser: haplotype‑aware SV genotyping},
  version      = {0.2.0},
  date         = {2024-06-18},
  url          = {https://github.com/your‑org/SvPhaser}
}
```




## License
`SvPhaser` is released under the MIT License – see [`LICENSE`](LICENSE).





## 📬 Contact

Developed by **Team5** (*BioAI Hackathon*) – Sachin Gadakh & Pranjul Mishra.

Lead contacts:
• [pranjul.mishra@proton.me](mailto:pranjul.mishra@proton.me)
• [s.gadakh@cent.uw.edu.pl](mailto:s.gadakh@cent.uw.edu.pl)

Feedback, feature requests and bug reports are all appreciated — feel free to open a GitHub issue or reach out by e‑mail.

---

*Happy phasing!*
