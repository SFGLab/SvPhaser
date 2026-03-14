# Methodology & Exploratory Results

This document provides a **rigorous, implementation‑faithful description** of the SvPhaser workflow *as implemented in v2.1.x*, together with guidance on how downstream benchmarking and interpretation should be performed. It reflects the **current architecture, decision logic, and output semantics** of SvPhaser after recent correctness fixes to GT / INFO serialization.

---

## 1  Pipeline overview

SvPhaser performs **post‑hoc haplotype assignment of pre‑called structural variants (SVs)** using phased long‑read alignments (HP tags in BAM/CRAM). The design goal is to remain **caller‑agnostic** while providing biologically interpretable haplotype labels and confidence metrics.

### 1.1 High‑level stages

1. **Input ingestion**

   * Unphased SV VCF (e.g. Sniffles2 output; optional `RNAMES` INFO field).
   * Phased long‑read BAM/CRAM with `HP` tags.

2. **Per‑chromosome evidence collection (parallelised)**

   * Each chromosome is processed independently by a worker.
   * For every SV, *ALT‑supporting reads* are identified using **SV‑type‑aware logic**:

     * **DEL**: large CIGAR `D` operations spanning breakpoints and/or split‑reads.
     * **INS**: large CIGAR `I` or soft‑clips near POS.
     * **BND**: split‑reads with SA tags linking partner breakpoints.
     * **INV**: split‑reads with strand inversion at breakpoints.

3. **Haplotype evidence aggregation**

   * Supporting reads are counted as:

     * `hp1`: reads with `HP=1`
     * `hp2`: reads with `HP=2`
     * `nohp`: supporting reads lacking HP tags
   * Derived quantities:

     * `tagged_total = hp1 + hp2`
     * `support_total = hp1 + hp2 + nohp`

4. **Haplotype classification (decision logic)**

   * A deterministic classifier assigns a phased genotype (`GT`) and genotype quality (`GQ`) using only **ALT‑support evidence**.
   * Classification is fully local (per‑SV) and independent across variants.

5. **Global filtering and serialization**

   * Variants with insufficient total support are dropped *after* per‑chromosome processing.
   * Final outputs:

     * `<stem>_phased.csv` (analysis‑first, row‑wise evidence table)
     * `<stem>_phased.vcf(.gz)` (VCF‑compliant GT + INFO annotations)

---

## 2  Decision model for haplotype assignment

SvPhaser deliberately separates **evidence collection** from **decision logic**. Only four quantities enter the classifier:

* `hp1`, `hp2`
* `support_total`
* user‑defined thresholds

### 2.1 Primary thresholds

| Parameter              | CLI flag               | Default | Interpretation                                                             |                   |             |
| ---------------------- | ---------------------- | ------- | -------------------------------------------------------------------------- | ----------------- | ----------- |
| Minimum support        | `--min-support`        | 10      | Drop SVs with insufficient total ALT support (HP + NOHP).                  |                   |             |
| Minimum tagged support | `--min-tagged-support` | 3       | Require at least this many HP‑tagged reads to attempt directional phasing. |                   |             |
| Major‑delta            | `--major-delta`        | 0.60    | Strong majority threshold: `max(hp1,hp2)/(hp1+hp2)`                        |                   |             |
| Equal‑delta            | `--equal-delta`        | 0.10    | Near‑tie threshold: `                                                      | hp1−hp2           | /(hp1+hp2)` |
| Tie handling           | `--tie-to-hom-alt`     | enabled | Near‑ties emitted as `1                                                    | 1`instead of`./.` |             |

### 2.2 Decision outcomes

Given sufficient support:

* **Directional heterozygous** (`1|0` or `0|1`)

  * `max(hp1,hp2)/(hp1+hp2) ≥ major_delta`

* **Homozygous alternate** (`1|1`)

  * `|hp1−hp2|/(hp1+hp2) ≤ equal_delta` and tie‑to‑hom enabled

* **Ambiguous** (`./.`)

  * Insufficient tagged reads or intermediate imbalance

Each decision is accompanied by:

* **GQ** – Phred‑scaled confidence derived from a binomial model on `(hp1, hp2)`.
* **REASON code** – categorical explanation (e.g. `MAJOR_HP1`, `TIE_HOM`, `LOW_TAGGED`).

---

## 3  Output semantics

### 3.1 CSV output (`*_phased.csv`)

The CSV is the **authoritative explainer** of SvPhaser behaviour. Each row corresponds to one retained SV and includes:

* Evidence counts: `hp1`, `hp2`, `nohp`, `tagged_total`, `support_total`
* Decision outputs: `gt`, `gq`, `reason`, `delta`
* Provenance: `mode`, `fetch_w`, `bp_tol`, `rnames_total`, `rnames_found`, `in_gt`

This table is intended for:

* threshold tuning
* statistical summaries
* debugging unexpected calls

### 3.2 VCF output (`*_phased.vcf`)

The phased VCF is **fully compliant** and contains:

* `FORMAT/GT` – phased genotype (`|` when applicable)
* `FORMAT/GQ` – integer genotype quality

If `--svp-info` is enabled, additional INFO fields are injected:

* `SVP_HP1`, `SVP_HP2`, `SVP_NOHP`
* `SVP_SUPPORT`, `SVP_TAGGED`, `SVP_TAGFRAC`
* `SVP_DELTA`, `SVP_REASON`
* `SVP_MODE`, `SVP_FETCHW`, `SVP_BPWIN`
* `SVP_RNAMES_TOTAL`, `SVP_RNAMES_FOUND`

Optionally, `GQBIN` labels are added if `--gq-bins` is provided.

---

## 4  Pseudoalgorithm (implementation‑faithful)

The following pseudoalgorithm mirrors the **actual control flow of SvPhaser v2.1.x**, abstracted from the implementation but preserving all decision points, thresholds, and invariants.

### 4.1 Per‑chromosome worker logic

```
for each chromosome C in genome (parallel):
    for each structural variant SV on chromosome C:

        # Step 1: identify ALT‑supporting reads
        if SUPPORT_MODE == RNAMES and SV has RNAMES:
            reads = fetch_reads_by_name(BAM, SV.RNAMES)
        else if SUPPORT_MODE == HEURISTIC:
            reads = fetch_reads_by_position(BAM, SV, window)
        else:  # HYBRID
            reads = RNAMES‑based if available else heuristic

        reads = filter_to_ALT_support(reads, SV.SVTYPE)

        # Step 2: aggregate haplotype evidence
        hp1  = count(read.HP == 1 for read in reads)
        hp2  = count(read.HP == 2 for read in reads)
        nohp = count(read.HP missing for read in reads)

        tagged_total  = hp1 + hp2
        support_total = hp1 + hp2 + nohp

        # Step 3: global support filter
        if support_total < MIN_SUPPORT:
            mark SV as DROPPED (LOW_SUPPORT)
            continue

        # Step 4: haplotype classification
        if tagged_total < MIN_TAGGED_SUPPORT:
            GT     = "./."
            GQ     = 0
            REASON = LOW_TAGGED
            DELTA  = 0.0

        else:
            major = max(hp1, hp2)
            minor = min(hp1, hp2)

            frac_major = major / tagged_total
            frac_diff  = (major - minor) / tagged_total

            if frac_major ≥ MAJOR_DELTA:
                GT     = "1|0" if hp1 > hp2 else "0|1"
                REASON = MAJOR_HP1 or MAJOR_HP2

            else if frac_diff ≤ EQUAL_DELTA:
                if TIE_TO_HOM_ALT:
                    GT     = "1|1"
                    REASON = TIE_HOM
                else:
                    GT     = "./."
                    REASON = TIE_AMBIG

            else:
                GT     = "./."
                REASON = AMBIG

            # Step 5: genotype quality
            GQ    = phred_scaled_binomial_confidence(hp1, hp2)
            DELTA = frac_diff

        # Step 6: emit row for CSV + VCF serialization
        emit_row(
            chrom, pos, id, svtype,
            hp1, hp2, nohp,
            tagged_total, support_total,
            GT, GQ, DELTA, REASON,
            fetch_window, bp_window,
            rnames_total, rnames_found,
            in_gt)
```

### 4.2 Global merge and serialization

```
merge rows from all chromosomes

assert column 'gt' exists and is populated

write CSV (analysis‑first artifact)

write VCF:
    FORMAT/GT = GT
    FORMAT/GQ = GQ
    INFO/SVP_* = evidence + provenance fields
```

This structure guarantees that **every phased VCF record is backed by an explicit CSV row**, and that no genotype is emitted without recorded evidence.

---

## 5  Workflow diagram

The overall SvPhaser workflow can be summarized as follows:

```
Unphased SV VCF        Phased BAM/CRAM
        │                     │
        └──────────┬──────────┘
                   ▼
        Per‑chromosome workers (parallel)
                   │
        ALT‑supporting read detection
                   │
        hp1 / hp2 / nohp aggregation
                   │
        Deterministic haplotype classifier
                   │
        (GT, GQ, REASON, DELTA)
                   │
        ┌──────────┴──────────┐
        ▼                     ▼
  Phased CSV           Phased VCF
 (analysis‑first)   (GT:GQ + SVP_*)
```

---

## 6  Interpretation guidance

* **CSV before VCF**: all benchmarking and biological reasoning should be driven from the CSV; the VCF is a consumable artifact.
* **1|1 calls** should be interpreted carefully and stratified by `reason` (true balance vs. tie‑forced homozygosity).
* **INS variants** are expected to show higher ambiguity rates due to alignment uncertainty.
* **NOHP fraction** is a key diagnostic for library or phasing quality issues.

---

## 5  Reproducibility checklist

* **Software:** SvPhaser v2.1.x
* **Execution model:** per‑chromosome parallel workers
* **Key inputs:** unphased SV VCF + HP‑tagged BAM
* **Command template:**

```bash
svphaser phase <sv.vcf> <phased.bam> \
  --out-dir results/ \
  --min-support 10 \
  --min-tagged-support 3 \
  --major-delta 0.60 \
  --equal-delta 0.10 \
  --support-mode hybrid \
  --dynamic-window \
  --tie-to-hom-alt \
  --gq-bins "30:High,10:Moderate" \
  --threads <N>
```

* **Reference genome:** same build used for alignment and SV calling.

---

*Last updated to reflect SvPhaser v2.1.x architecture and serialization guarantees.*
