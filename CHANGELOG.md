# Changelog

All notable changes to **SvPhaser** will be documented in this file. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project adheres to [Semantic Versioning](https://semver.org/).

---

## \[0.2.0] – 2025‑06‑18

### Added

* **Parallel worker layer** (`svphaser/phasing/_workers.py`) that cleanly isolates per‑chromosome work and is pickle‑safe for `multiprocessing`.
* **Typed dataclasses & aliases** in `svphaser/phasing/types.py` for stronger static typing across the code‑base.
* **Logging helper** (`svphaser/logging.py`) with coloured, timestamped output and per‑module log‑level control.
* **GQ binning** via new `--gq-bins` CLI option and associated parsing/annotation logic.
* Rich documentation: `docs/methodology.md` now contains an updated Mermaid workflow diagram plus exploratory result plots.
* Example result plots & workflow PNGs under `docs/result_images/` (CI ignores via `.gitattributes`).

### Changed

* Refactored `svphaser/phasing/io.py` to:

  * Use `multiprocessing.Pool.starmap` for clearer argument passing.
  * Replace cyvcf2 `VCF.copy()` (deprecated) with manual header duplication to maintain compatibility.
  * Centralise depth filtering logic and move decision thresholds into `WorkerOpts`.
* CLI (`svphaser/cli.py`):

  * Switched `--out-vcf` & `--out-csv` to a single **`--out-dir`** argument that auto‑generates stemmed filenames.
  * Added graceful error trapping with coloured Typer feedback.

### Fixed

* Type‑checking errors flagged by Pylance/MyPy (e.g. `Dict[SVKey, …]` mismatches, `bam_path: Path` incompatibility, missing `gq_bins`).
* Multiprocessing pickling crash caused by passing `cyvcf2.Variant` directly.
* AttributeError for `cyvcf2.VCF.copy()` on older cyvcf2 builds – now handled by explicit header cloning.
* Jupyter export issue for SVG → PNG in VS Code notebooks (`Unsupported MimeType`) by switching to `plt.savefig(..., format="png")`.

### Removed

* Obsolete helper functions in `algorithms.py` superseded by the new decision‑tree implementation.

---

## \[0.1.0] – 2025‑05‑15

### Added

* First public release with core phasing algorithm, CSV/VCF output and basic unit tests.
