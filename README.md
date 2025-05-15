# SvPhaser

[![PyPI version](https://badge.fury.io/py/svphaser.svg)](https://badge.fury.io/py/svphaser)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Structural Variant Phasing Tool for Long-reads based Phased BAM - Developed by Team5 (BioAI Hackathon)

---

## 🧬 Overview

**SvPhaser** is a lightweight, fast, and scalable tool for **phasing structural variants (SVs)** in human genomes by leveraging:
- **Phased aligned BAM files** (containing HP tags)
- **Unphased structural variant VCF files**

Unlike traditional SV callers, **SvPhaser does not call new SVs** — it assigns haplotypes to **existing SVs** based on read evidence.

---

## 📈 Pipeline Workflow

![Phasing SV Pipeline](Output/Pipeline_Diagram.png)

```
Input: Phased BAM + Unphased SV VCF
↓
Parallel chromosome-wise read extraction
↓
Read support evaluation & haplotype assignment
↓
Merged CSV creation
↓
Phased SV VCF writing
↓
Final Output
```
## 📝 Our Approach

- Designed to assign haplotypes to pre-detected structural variants (SVs) using phased long-read BAM alignments.
- Extracts HP-tagged reads overlapping SV regions.
- Assigns SVs to haplotype 1 or 2 by majority of supporting reads.
- Marks SV as HP1_HP2 (CSV) / 1|1 (VCF) when equal support exists.
- Marks SV as unphased (unphased / ./. in VCF) if below read support threshold.
- Runs efficiently in parallel on chromosome-split datasets.
- Outputs per-chromosome CSVs, a merged CSV, and final phased VCF file.

---

## 🚀 Key Features

- Parallel chromosome-wise SV phasing
- Flexible minimum read support threshold
- Outputs phased SV VCF with GT and DP fields
- Fast processing with low memory usage
- Pip-installable + easy CLI interface

---

## ⚙️ System Requirements

- OS: Linux / WSL2 (Ubuntu 22.04 recommended)
- Python: >= 3.6
- CPU: 4 cores minimum (parallelization supported)
- RAM: 16 GB minimum (recommended 24 GB+ for full genomes)

---

## 🛠️ Installation

### 1. Clone repository
```bash
git clone https://github.com/SFGLab/SvPhaser.git
cd SvPhaser/SvPhaser
```

### 2. Install via pip
```bash
pip install .
```

or for development mode:
```bash
pip install -e .
```

---

## 📦 Dependencies

| Package | Purpose |
|---------|---------|
| pysam | BAM file parsing |
| pandas | VCF/CSV parsing |
| argparse | CLI parsing (built-in) |
| glob, multiprocessing, os, collections | Built-in libraries |

✅ Only `pysam` and `pandas` are external.

---

## 🚀 Usage Example

```bash
svphaser --phased_bam HG00733.sorted_phased.bam \
         --unphased_vcf HG00733_allsvs_10X.vcf \
         --output /path/to/output_folder \
         --min_support 10
```

**Arguments:**

| Parameter | Description |
|-----------|-------------|
| `--phased_bam` | Path to BAM file containing HP tags |
| `--unphased_vcf` | Path to unphased structural variant VCF |
| `--output` | Output directory |
| `--min_support` | Minimum number of supporting reads (default: 10) |

---

## 📄 Output Structure

```
output_folder/
├── chromosome_csvs/ (temporary CSVs - deleted automatically)
├── merged/
│   ├── SV_phasing_full.csv
│   ├── {input_vcf_name}_phased.vcf
```


## 📊 Results of Our Analysis

### 1️⃣ Overlap of Phased SV Sets

This pie chart summarizes overlap between diploid-phased SVs and alignment-based phased SVs using a 1bp overlap threshold.

![Overlap of Phased SV Sets](Output/graphs/overlap_phased_sv_sets_piechart.png)

### 2️⃣ IGV Validation of Phased SV

An example IGV screenshot showing an unphased SV (original), the same SV successfully phased by SvPhaser, and confirmation in the diploid assembly.

![IGV Validation of Phased SV](Output/graphs/igv_phased_sv_validation_chr16.png)

### 3️⃣ Haplotype Assignment Distribution

Bar plot showing SV assignment frequencies by haplotype with min_support=10. Majority of SVs are confidently assigned to HP1 or HP2.

![SV Phasing Results](Output/graphs/sv_phasing_results_min_support_10_barplot.png)

![Distribution Description](Output/graphs/Distribution_Description.png) 


### How will VCF Phased Output look like ? 
We made some ammenmends and made the vcf file more informative which is compatible with IGV so no error with the changes 
![Phased VCF Output File](Output/graphs/Vcf_output_phased.png) 
 
---

## 📊 Benchmarking

| System | CPU | RAM | Runtime | Dataset Size |
|--------|-----|-----|---------|--------------|
| Workstation (Linux) | 14 cores | 256 GB | ~30 minutes | Full human genome (30x) |
| Laptop (WSL2, i7-9th Gen) | 6 cores | 24 GB | ~1hours | Full human genome (30x) |

✅ Scales linearly with number of CPU cores.  
✅ Fully memory-safe even on small systems.

---


## 💎 Novel Contributions

- Direct haplotype assignment to SVs using phased BAM read evidence
- Per-SV read support statistics integrated in VCF output
- Automatic deletion of intermediate files for disk efficiency
- Lightweight CLI designed for both local and cluster use
- Ready for integration into large cohort SV analysis pipelines

---

## 📜 License

This project is licensed under the [MIT License](LICENSE).

---

## 📬 Contact

Developed by: Team5 (BioAI Hackathon): Sachin Gadakh, Pranjul Mishra  
Lead Contacts: [pranjul.mishra@proton.me], [s.gadakh@cent.uw.edu.pl]  

Feel free to submit issues, feature requests, or contribute via GitHub!

---

# 🚀 Ready to Phase Your SVs?
```bash
pip install svphaser
svphaser --help
```

---

# 🎉 Release Notes (v1.0.0)

- Initial public release of **SvPhaser**
- Full parallel chromosome-wise phasing
- CLI interface with `--min_support` option
- Integrated VCF writer with read support field
- Disk-efficient: deletes temporary files after merge

---
