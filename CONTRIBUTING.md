# Contributing to **SvPhaser**

Thanks for your interest in improving SvPhaser! 🎉
This guide explains how to set up a development environment, our coding standards, test/benchmark expectations, and how to submit changes.

> By participating, you agree to follow our [Code of Conduct](CODE_OF_CONDUCT.md).

---

## Table of contents

* [Quick start](#quick-start)
* [Project layout](#project-layout)
* [Development workflow](#development-workflow)
* [Coding standards](#coding-standards)
* [Testing](#testing)
* [Performance & profiling](#performance--profiling)
* [VCF & I/O validation](#vcf--io-validation)
* [Documentation](#documentation)
* [Git & PR guidelines](#git--pr-guidelines)
* [Release process](#release-process)
* [Issue reports & feature requests](#issue-reports--feature-requests)
* [Security](#security)
* [License](#license)

---

## Quick start

**Prerequisites**

* Python **3.10+**
* `pip`, `virtualenv` (or `conda`/`mamba`)
* Development tools: a C/C++ compiler (for `pysam`/`cyvcf2`), `htslib` headers recommended

```bash
# 1) Create and activate a virtualenv
python -m venv .venv
source .venv/bin/activate     # Windows: .venv\Scripts\activate

# 2) Upgrade packaging tools
python -m pip install -U pip wheel

# 3) Install SvPhaser in editable mode (runtime deps from pyproject/requirements)
pip install -e .

# 4) Install development deps
pip install -r requirements-dev.txt

# 5) (Recommended) Install pre-commit hooks
pre-commit install

# 6) Run tests to verify your setup
pytest -q
```

**Running the CLI locally**

```bash
svphaser <SV_VCF> <HP-tagged BAM/CRAM> \
  --out-dir out/ --min-support 10 --major-delta 0.60 --equal-delta 0.10 \
  --gq-bins "30:High,10:Moderate" --threads 8
```

> **Large files**: please do **not** commit BAM/CRAM/large VCFs. For tests, use tiny synthetic fixtures (≤200 KB) under `tests/data/`. For integration benchmarks, point tests to external paths via env vars (see below).

---

## Project layout

```
src/svphaser/
  __init__.py          # public API: phase() wrapper, version, defaults
  __main__.py          # entry point for python -m svphaser
  cli.py               # Typer CLI app (phase subcommand)
  logging.py           # module-level logging configuration
  phasing/
    __init__.py        # exports: phase_vcf, classify_haplotype, phasing_gq, WorkerOpts
    algorithms.py      # pure math: phasing_gq(), classify_haplotype() (no I/O)
    io.py              # orchestration: VCF parsing, worker spawning, CSV/VCF writing
    _workers.py        # internal: _phase_chrom_worker() (per-chromosome logic, BAM parsing)
    types.py           # dataclasses: WorkerOpts, NamedTuple: CallTuple; type aliases

tests/
  test_algorithms.py   # GQ monotonicity, symmetry, symmetry; near-threshold edge cases
  test_cli_smoke.py    # CLI parsing, help text
  test_io.py          # CSV column order, VCF header preservation
  test_workers.py     # HP tag filtering, size consistency, read counting
```

### Code organization principles

* **Pure code vs. I/O**: Keep probability/math in `algorithms.py` (no file I/O). Keep BAM/VCF handling in workers/IO.
* **Parallelism**: Multiprocessing by chromosome; all inputs must be pickle-friendly.
* **Determinism**: No random sampling; set seeds explicitly in hypothesis tests.
* **Typing**: All public functions must have type hints; `mypy --strict` on `src/`.
* **Logging**: Use module loggers (e.g., `logging.getLogger(__name__)`); no prints in library code.

---

## Development workflow

1. **Create a branch**

   * `git switch -c feature/<topic>` or `fix/<bug>`
2. **Write the change**

   * Add/adjust unit tests first when you can (TDD encouraged).
3. **Run quality checks**

   ```bash
   ruff check src tests
   black --check src tests
   mypy src
   pytest -q
   ```
4. **Commit using Conventional Commits** (see below) and push.
5. **Open a Pull Request** (PR) with a clear description and checklist.

### Conventional Commits

Examples: `feat: add --regions filter`, `fix(io): tab-delimit VCF lines`, `perf(workers): reduce BAM scans`.

---

## Coding standards

* **Types**: use **type hints** on all public functions; use generics where applicable.
* **Formatting**: `black` (line length 100, Python 3.9 target).
* **Linting**: `ruff check` for E/F/W/I/UP/B/C90; fix or justify inline.
* **Static typing**: `mypy --strict` on `src/` (pandas-stubs + hypothesis stubs allowed).
* **Docstrings**: concise, with argument/return docs on non-trivial functions (NumPy style).
* **Logging**: module loggers only (`logging.getLogger(__name__)`); no `print()` statements.
* **Error handling**: raise informative `ValueError`/`RuntimeError` with actionable messages.
  - Bad: `raise ValueError("invalid input")`
  - Good: `raise ValueError(f"BAM has no index; please run 'samtools index {bam_path}'")`
* **Determinism**: tests must be deterministic; use `@given` seeds and track randomness.
* **Module constants**: centralize default thresholds in `__init__.py` (e.g., `DEFAULT_MIN_SUPPORT = 10`).
* **Comments**: explain *why*, not *what*; complex logic gets a 1-2 line docstring.

### Example structure for a new function

```python
def estimate_sv_coverage(
    bam_path: Path,
    region: str,
    *,
    sample_rate: float = 0.1,
    min_mapq: int = 20,
) -> int:
    """Estimate mean SV-supporting read depth in a region.

    Args:
        bam_path: Path to coordinate-sorted, indexed BAM/CRAM.
        region: Region string (e.g., "chr1:1000-2000").
        sample_rate: Fraction of reads to sample [0.0, 1.0].
        min_mapq: Skip reads with MAPQ < this.

    Returns:
        Mean number of supporting reads per position.

    Raises:
        FileNotFoundError: if BAM or index is missing.
        ValueError: if sample_rate outside [0, 1].
    """
    if not 0.0 <= sample_rate <= 1.0:
        raise ValueError(f"sample_rate must be in [0.0, 1.0]; got {sample_rate}")
    # ... implementation
```

**PR checklist (copy into PR body):**

* [ ] Unit tests added/updated for new logic
* [ ] Lint/format pass (`ruff check src tests`, `black src tests`)
* [ ] Type-check pass (`mypy src`)
* [ ] VCF validation pass (see § VCF & I/O validation)
* [ ] Benchmarks attached or old numbers confirmed (if perf-sensitive)
* [ ] Docs/CLI help updated (README, docstrings, `--help` output)
* [ ] CHANGELOG.md entry added (for user-facing changes)

---

## Testing

We use **pytest**. Small test data lives in `tests/data/` (minified VCF/BAM slices). For large integration tests use env variables:

* `SVPHASER_TEST_VCF` — path to a real VCF (optional)
* `SVPHASER_TEST_BAM` — path to a real HP-tagged BAM/CRAM (optional)

```bash
pytest -q                              # quick run
pytest -q -n auto                      # parallel (requires pytest-xdist)
pytest --cov=svphaser --cov-report=term-missing  # coverage
pytest tests/test_algorithms.py -v     # verbose
pytest tests/test_io.py --hypothesis-seed=12345  # reproducible Hypothesis runs
```

### What to test

**`algorithms.py` (core math)**
* GQ monotonicity: `GQ(n1, n2)` ≥ `GQ(n1-1, n2)` for same N
* GQ symmetry: `GQ(a, b) == GQ(b, a)` and same for delta
* GQ capping: `GQ(...) ≤ 99` always
* Edge cases: `n1=0, n2=0`, `n1>>n2`, deep coverage (N>200, normal approx)
* Classification thresholds: test near `major_delta`, `equal_delta`, `min_support` boundaries
* Tie handling: verify `tie_to_hom_alt=True` gives `1|1`, `False` gives `./.`

**`_workers.py` (BAM parsing)**
* HP tag filtering: correctly identifies `HP=1`, `HP=2`, missing HP
* Secondary/supplementary filtering: unmapped/secondary/supplementary reads excluded
* Size consistency (DEL/INS): variant size matches computed from reads ± tolerance
* Multi-allelic handling: correct ALT indexing for multiple alts
* Reason codes: "MinSupport", "LowTagged", "Tie", "Size-mismatch" etc.

**`io.py` (CSV/VCF writing)**
* CSV headers present: `chrom`, `pos`, `id`, `end`, `svtype`, `gt`, `gq`, etc.
* CSV tab-delimited: no extra spaces
* VCF header preservation: original `#CHROM` line unchanged
* VCF `FORMAT=GT:GQ`: correct sample columns and value formats
* Optional columns: `gq_label` computed correctly if `gq_bins` provided
* `tag_frac` backfill: `tagged_total / support_total` (or NaN if zeros)
* Dropped SVs: written to `*_dropped_svs.csv` with reason

**`cli.py` (command-line)**
* All required args parsed: VCF, BAM, out-dir
* Optional flags work: `--threads`, `--gq-bins`, `--no-svp-info`, etc.
* Help text accurate: `svphaser --help`, `svphaser phase --help`
* Exit codes: 0 on success, ≠0 on error
* Error messages helpful: "File not found: ...", "Invalid --support-mode ..."

**Example test patterns**

```python
from hypothesis import given, strategies as st
import pytest

def test_gq_monotonicity():
    """GQ increases with |n1-n2| at fixed N."""
    gq1 = phasing_gq(8, 2)
    gq2 = phasing_gq(9, 1)
    assert gq2 >= gq1

@given(st.integers(0, 500), st.integers(0, 500))
def test_gq_symmetry(n1, n2):
    """GQ(a, b) == GQ(b, a)."""
    assert phasing_gq(n1, n2) == phasing_gq(n2, n1)

def test_classify_tie_with_tags():
    """When |n1-n2| ≤ equal_delta and tie_to_hom_alt=True, emit 1|1."""
    gt, gq = classify_haplotype(
        5, 4, equal_delta=0.2, tie_to_hom_alt=True, min_support=9
    )
    assert gt == "1|1"

def test_cli_version():
    """svphaser --version prints version and exits."""
    from svphaser.cli import app
    from typer.testing import CliRunner
    runner = CliRunner()
    result = runner.invoke(app, ["--version"])
    assert result.exit_code == 0
```

---

## Performance & profiling

Before merging changes that touch hot paths (BAM parsing, GQ calculation, worker orchestration):

```bash
# Quick benchmark
time svphaser phase sample.vcf.gz sample.bam --threads 1 -o out/
time svphaser phase sample.vcf.gz sample.bam --threads 8 -o out/

# Memory profile
/usr/bin/time -v svphaser phase sample.vcf.gz sample.bam -o out/
# Look for: "Maximum resident set size"

# CPU profile (requires py-spy, install with pip install py-spy)
py-spy record -o prof.svg -- svphaser phase sample.vcf.gz sample.bam -o out/

# Python cProfile
python -m cProfile -s cumulative -m svphaser phase sample.vcf.gz sample.bam -o out/ \
  | head -30

# Check scaling with threads
for n in 1 2 4 8 16; do
  echo "threads=$n:"
  time svphaser phase sample.vcf.gz sample.bam --threads $n -o out_$n/ 2>&1 | grep real
done
```

**What to watch for**

* **BAM I/O**: ensure tabix-indexed VCF allows region iteration (not full scan).
  - Check: `cyvcf2.Reader` is opened with `region` parameter (see `io.py`).
* **GQ calculation**: binomial tail up to N=200, then Gaussian. Ensure no overflow at huge depths.
* **Worker pickling**: multiprocessing workers must serialize inputs; avoid large unpickleable objects.
  - Fail early: `WorkerOpts` is a frozen dataclass ✓; don't pass open file handles to workers ✗.
* **Memory**: RSS should stay ~constant across chromosomes (not accumulate). Watch for pandas `concat()` leaks.
* **Parallel speedup**: expect near-linear speedup for N threads up to # CPUs. Beyond, contention rises.

Document any regressions/improvements in the PR (e.g., "Reduced peak RSS from 2.1 GB → 1.8 GB on chr22", "Threadpool now scales 7.8× on 8 cores").

---

## VCF & I/O validation

Every PR that touches VCF/CSV/BAM logic should validate outputs:

```bash
# Parse & check headers
bcftools view -h out/sample_phased.vcf | head -5
bcftools view -h out/sample_phased.vcf | grep "^#CHROM"  # sample names present

# Spot-check records
bcftools view out/sample_phased.vcf | head -3

# Validate VCF compliance (if vcftools/vt available)
bcftools view out/sample_phased.vcf > /tmp/test.vcf  # decompress to plain VCF
vcf-validator /tmp/test.vcf

# Check CSV integrity
head tests/data/sample_phased.csv  # headers present
tail tests/data/sample_phased.csv  # no truncation
wc -l tests/data/sample_phased.csv  # same as VCF variant count (±1 for header)

# Spot-check data alignment
cat tests/data/sample_phased.csv | cut -f1-3,6,7  # chrom, pos, id, gt, gq
```

### VCF conventions we follow

* **Headers**: preserve original `#CHROM` line; add SvPhaser version & command line as comment.
* **Tab-delimited**: all fields separated by single `\t` (no spaces around delimiters).
* **INFO preservation**: pass through original `REF/ALT/QUAL/FILTER/INFO` unchanged.
  - Prepend `SVTYPE` (from ALT tag or infer) early in INFO.
  - Append `SVP_HP1`, `SVP_HP2`, `SVP_NOHP`, `SVP_TAGFRAC`, `SVP_DELTA`, `SVP_GQBIN` (if requested).
* **FORMAT**: always `GT:GQ` (genotype, then quality).
  - `GT` values: `0/1`, `1/1`, `./.`, never `0/0` (SV neutral genotype is `./.` not `0/0`).
  - `GQ` values: integer 0–99.
* **Sample column**: derive sample name from input VCF header (or if missing, "sample").

### CSV conventions we follow

* **Tab-delimited** (like VCF).
* **Must include**: `chrom`, `pos`, `id`, `end`, `svtype`, `gt`, `gq`.
* **Optional but common**: `hp1`, `hp2`, `nohp`, `tagged_total`, `support_total`, `delta`, `equal_delta`, `tag_frac`, `gq_label`, `reason`.
* **No sorting required** (output order follows input VCF + workers).
* **Dropped SVs**: written to separate `*_dropped_svs.csv` with reason codes.

### Example validation in tests

```python
import pandas as pd
from pathlib import Path

def test_output_csv_format():
    """CSV is tab-delimited with required columns."""
    csv_path = Path("results/sample_phased.csv")
    df = pd.read_csv(csv_path, sep="\t")

    required = {"chrom", "pos", "id", "end", "svtype", "gt", "gq"}
    assert required.issubset(df.columns), f"Missing: {required - set(df.columns)}"

    # Validate GT values
    assert df["gt"].isin(["0/1", "1/0", "1/1", "./."]).all(), "Invalid GT values"

    # Validate GQ range
    assert (df["gq"] >= 0).all() and (df["gq"] <= 99).all(), "GQ out of range [0, 99]"
```

---

## Documentation

* Keep **README.md** and CLI `--help` text in sync (update both when adding flags).
* **Algorithm changes**: update `docs/Methodology.md` and/or inline docstrings with rationale.
* **Output format changes**: document in README outputs section and CSV/VCF columns.
* **Large changes**: consider adding a subsection in docs/ (e.g., `docs/SIZE_MATCHING.md`).
* **Code comments**: explain *why*, not *what*. Good: `# Emit 1|1 when tied support but both haplotypes present`; Bad: `# Set gt to 1|1`.
* **Docstrings**: use NumPy style, include examples for public functions.

Example docstring:

```python
def classify_haplotype(...) -> tuple[str, int]:
    """Assign haplotype-aware genotype to an SV.

    Semantics follow SvPhaser v2.1.1 decision tree:
    1) If total support < min_support: emit ./.
    2) If tagged < min_tagged_support: emit ./.
    3) If |HP1−HP2| / tagged ≤ equal_delta: emit 1|1 (if both > 0 and tie_to_hom_alt) else ./.
    4) Else if max(HP1, HP2) / tagged ≥ major_delta: emit 1|0 or 0|1 (directional).
    5) Else: emit ./.

    Parameters
    ----------
    n1 : int
        Number of reads tagged HP=1.
    n2 : int
        Number of reads tagged HP=2.
    min_support : int, default 10
        Minimum total ALT-supporting reads (tagged + untagged).
    ...

    Returns
    -------
    tuple[str, int]
        (genotype, gq) where genotype ∈ {"1|0", "0|1", "1|1", "./."} and
        gq is Phred-scaled quality (0–99).

    Examples
    --------
    >>> classify_haplotype(20, 5, min_support=10, major_delta=0.6)
    ('1|0', 45)
    """
```

---

## Git & PR guidelines

### Branch naming

* `feature/<topic>` — new capability (e.g., `feature/regions-filter`)
* `fix/<bug>` — bug fix (e.g., `fix/haplotype-tie-logic`)
* `perf/<area>` — performance improvement (e.g., `perf/bam-io-caching`)
* `refactor/<area>` — code reorganization without behavior change (e.g., `refactor/algorithm-docs`)
* `test/<topic>` — test additions/improvements (e.g., `test/edge-case-coverage`)

### Commit messages

Use [Conventional Commits](https://www.conventionalcommits.org/):

```
feat: add --regions filter for targeting specific genomic areas
fix(algorithms): correct GQ calculation for near-zero tail probability
perf(workers): optimize BAM seek overhead with tabix indexing
docs: expand README with parameter table
test: add property-based tests for GQ symmetry
```

Not: "updated code", "fixes", "wip".

### Pull request checklist

Before opening a PR, ensure:

- [ ] **Branch is up-to-date** with `main` / `develop`
- [ ] **Tests pass**: `pytest -q`
- [ ] **Linting passes**: `ruff check src tests && black --check src tests`
- [ ] **Type checking passes**: `mypy src`
- [ ] **VCF/CSV validation** runs (see § VCF & I/O validation)
- [ ] **CHANGELOG.md** entry added (if user-facing)
- [ ] **Documentation** updated (README, docstrings, `--help` output)
- [ ] **Commit messages** follow Conventional Commits
- [ ] **PR description** explains *why* (linked issue, design rationale)

### PR description template

```markdown
## What

Brief summary of changes (1–2 sentences).

## Why

Motivation: what problem does this solve? Link to issue if applicable: Fixes #123.

## How

Implementation approach. Highlight any trade-offs or design decisions.

## Testing

How was this tested? If benchmarks affected, include before/after timings.

## Checklist

- [x] Tests added/updated
- [x] Lint/format pass
- [x] Type-check pass
- [x] VCF/CSV validation
- [x] Docs updated
```

### Review process

* At least **one approval** required before merge.
* Reviewers check: correctness, performance, test coverage, documentation, API compatibility.
* Address feedback or ask clarifying questions.
* Squash & merge on approval (one commit per feature).

---

## Release process

* We follow **semantic versioning** (MAJOR.MINOR.PATCH).
* Bump version in `pyproject.toml` (and `svphaser/phasing/__init__.py` if it exposes `__version__`).
* Update `CHANGELOG.md`.
* Build & upload:

```bash
python -m build
python -m twine upload dist/*
```

* Create a GitHub release with highlights and checksums.

---

## Issue reports & feature requests

### Reporting a bug

Include:

* **Version**: `svphaser --version`
* **OS**: `uname -a` or Windows version
* **Python version**: `python --version`
* **Exact command** you ran (sanitized, no credentials)
* **Full error message** and traceback
* **Minimal reproducible example**:
  - If file-based: share or link to tiny test data (≤1 MB VCF + BAM slice)
  - If numeric: provide exact `--min-support`, `--major-delta` etc.
* **Expected behavior** vs. **actual behavior**
* **Environment variables**: `SVPHASER_*` flags you set, if any

Example:

```
## Reproducibility

- SvPhaser version: 2.1.3
- Python 3.11.4 on Ubuntu 22.04
- Command:
  ```bash
  svphaser phase sample.vcf.gz sample.bam --threads 8 --major-delta 0.5
  ```
- Error:
  ```
  ValueError: BAM file has no index. Run: samtools index sample.bam
  ```

## Expected behavior

After indexing, should generate phased.vcf and phased.csv.

## Actual behavior

Crashes even after indexing (see attached samtools stderr).
```

### Requesting a feature

Describe:

* **Use case**: why would you need this? What workflow does it enable?
* **Expected API/CLI**: Show an example of how you'd use it.
* **Interaction with existing flags**: Does it conflict with `--tie-to-hom-alt` or other params?
* **Output impact**: Does it change CSV columns, VCF INFO tags, etc.?

Example:

```
## Feature: --min-allele-count filter

Currently, SvPhaser phases all SVs ≥ --min-support.
I need to filter down to variants present in ≥2 samples (i.e., n_samples ≥ 2).

Proposed CLI:
  --min-allele-count 2  (or --mac 2)

This would count how many samples carry the ALT at this locus and drop variants below threshold.

Output: add MAC column to CSV.
```

---

## Security

If you discover a security or privacy issue:

* **Do not** post it as a public GitHub issue.
* **Email** the maintainers privately: **[pranjul.mishra@proton.me](mailto:pranjul.mishra@proton.me)**
* Include: affected versions, detailed reproduction steps, potential impact.
* We will respond within 48 hours and coordinate a fix privately before public disclosure.

Common security concerns in bioinformatics:

* **PHI/PII in test data**: never commit real patient samples or identifiers
* **Temporary files**: ensure `.tmp.vcf` is deleted even on crash (use `tempfile` or context managers)
* **Command injection**: always quote filenames/paths; use `subprocess.run(..., shell=False)` or `pathlib.Path`
* **Input validation**: reject malformed VCF/BAM early with clear error messages

---

## License

By contributing, you agree that your contributions will be licensed under the project’s license (see `LICENSE`).
