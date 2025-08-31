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
  phasing/
    algorithms.py   # pure math: GQ, decision tree
    _workers.py     # per-chromosome worker, read counting
    io.py           # orchestration, CSV/VCF writers
    types.py        # shared light datatypes
    __init__.py     # public API exports
```

* **Pure code vs. I/O**: keep probability/math in `algorithms.py` (no file I/O); keep BAM/VCF handling in workers/IO.
* **Parallelism**: multiprocessing by chromosome (process-safe, pickle-friendly inputs only).

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

* **Types**: use **type hints** everywhere; keep public functions typed.
* **Formatting**: `black` (line length 88).
* **Linting**: `ruff` for style & imports; fix warnings or add brief justifications.
* **Static typing**: `mypy` on `src/` (allow `typing`/`pandas-stubs`).
* **Docstrings**: concise, with arguments/returns on non-trivial functions.
* **Logging**: module loggers only (no prints), e.g. `logging.getLogger("svphaser.io")`.
* **Error handling**: clear errors with actionable messages (e.g., missing index, bad header).
* **Determinism**: tests must be deterministic; set seeds where randomness exists.

**PR checklist (copy into PR):**

* [ ] Unit tests added/updated
* [ ] Lint/format pass (`ruff`, `black`)
* [ ] Type-check pass (`mypy`)
* [ ] VCF validation pass (see below)
* [ ] Benchmarks (if perf-sensitive) attached or unchanged
* [ ] Docs/CLI help updated

---

## Testing

We use **pytest**. Small test data lives in `tests/data/` (minified VCF/BAM slices). For large integration tests use env variables:

* `SVPHASER_TEST_VCF` — path to a real VCF (optional)
* `SVPHASER_TEST_BAM` — path to a real HP-tagged BAM/CRAM (optional)

```bash
pytest -q
pytest -q -n auto                      # parallel (if pytest-xdist installed)
pytest --cov=svphaser --cov-report=term-missing
```

**What to test**

* Decision tree edges: near thresholds (e.g., `major_delta`, `equal_delta`), low depth, deep coverage (GQ capping), ties.
* BAM read counting: ignores unmapped/secondary/supplementary; HP tag logic.
* CSV/VCF writing: headers, tab-delimited, sample column order.
* CLI: required/optional flags, helpful errors.

**Property tests** (optional): with Hypothesis, exercise GQ monotonicity, symmetry, and bounds.

---

## Performance & profiling

Before merging changes that affect hot paths:

* Time with `/usr/bin/time -v` and record RSS/CPU.
* Profile CPU: `python -m cProfile -o prof.out -m svphaser ...`; inspect with `snakeviz prof.out`.
* Track parallel speedup with `--threads` (1, 2, 4, …). Ensure no pickling of non-serializable objects.
* Watch I/O: ensure tabix-indexed VCF gets region iteration; fallback linear scan only when needed.

Document any regressions/improvements in the PR.

---

## VCF & I/O validation

Every PR that touches VCF/CSV/BAM logic should run at least:

```bash
# Header & parse check
bcftools view -h out/sample_phased.vcf > /dev/null
bcftools view out/sample_phased.vcf | head

# Optional deep validation if available
# vt validate out/sample_phased.vcf
# vcf-validator out/sample_phased.vcf
```

**VCF conventions we follow**

* Write the original header + a proper `#CHROM` header line.
* Tab-delimit all fields (no spaces).
* Preserve original `REF/ALT/QUAL/FILTER/INFO`; add `SVTYPE` (first in INFO), `GQBIN` if present.
* `FORMAT` = `GT:GQ`; sample name derived from input header; GQ capped at 99.

---

## Documentation

* Keep the README and CLI `--help` accurate.
* If you add flags/outputs, update examples and docstrings.
* Longer docs (design, benchmarks) live under `docs/` (if present) or as Markdown files at repo root.

---

## Git & PR guidelines

* Branch naming: `feature/<topic>`, `fix/<bug>`, `perf/<area>`, `refactor/<area>`.
* Keep PRs focused (< \~500 LOC when possible). Split refactors from behavior changes.
* Include **before/after** numbers for performance-affecting code.
* Mark **breaking changes** clearly and update the changelog.

**Reviews**

* At least one approval required.
* Reviewers check: correctness, performance risks, tests, docs, and user-facing compatibility.

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

When filing an issue, include:

* SvPhaser version (`svphaser --version`), OS, Python version
* Exact command you ran, full log (add `-v` for verbosity)
* Minimal input snippet if possible (VCF header + 1–3 records)
* Expected vs. actual behavior
* For performance issues: dataset size, `--threads`, `/usr/bin/time -v` output

Feature requests welcome—explain the use case, expected CLI/outputs, and how it interacts with existing flags.

---

## Security

If you discover a security or privacy issue, please contact the maintainers privately at **[pranjul.mishra@proton.me](mailto:add-contact@your.org)**. We will respond promptly and coordinate a fix.

---

## License

By contributing, you agree that your contributions will be licensed under the project’s license (see `LICENSE`).
