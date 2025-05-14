import pysam
import pandas as pd
from collections import defaultdict, Counter
import os

print("Libraries loaded ✅")

# Customize this based on where your drive is mounted
bam_file = "/media/your_drive_name/your_dataset_folder/HG00733.sorted_phased.bam"
sv_vcf = "/media/your_drive_name/your_dataset_folder/whatshap_phased.vcf.gz"
output_file = "/media/your_drive_name/your_dataset_folder/output_sv_phasing.csv"

print("Paths set ✅")

def get_read_haplotypes(bam_file):
    haplotypes = defaultdict(list)
    bam = pysam.AlignmentFile(bam_file, "rb")
    for read in bam.fetch():
        if read.has_tag("HP"):
            hp_tag = read.get_tag("HP")
            haplotypes[read.reference_name].append((read.reference_start, read.reference_end, hp_tag))
    bam.close()
    return haplotypes

def assign_sv_to_haplotype(sv_vcf, haplotypes, output_file, min_support=5):
    compression = 'gzip' if sv_vcf.endswith('.gz') else None

    sv_df = pd.read_csv(
        sv_vcf,
        comment="#",
        sep="\t",
        header=None,
        compression=compression
    )

    sv_df.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]

    results = []

    for _, row in sv_df.iterrows():
        chrom = row["CHROM"]
        pos = int(row["POS"])
        sv_start = pos
        sv_end = pos + 1

        haps = [
            hp for start, end, hp in haplotypes.get(chrom, [])
            if start <= sv_start <= end
        ]

        support_count = len(haps)
        assigned_hap = "unknown"

        if support_count >= min_support:
            hap_counter = Counter(haps)
            hap_freq = hap_counter.most_common()

            if len(hap_freq) >= 2 and hap_freq[0][1] == hap_freq[1][1]:
                assigned_hap = "HP1_HP2"
            else:
                if hap_freq[0][0] == 1:
                    assigned_hap = "HP1"
                elif hap_freq[0][0] == 2:
                    assigned_hap = "HP2"

        results.append([chrom, pos, row["ID"], assigned_hap, support_count])

    results_df = pd.DataFrame(results, columns=["CHROM", "POS", "SV_ID", "HAPLOTYPE", "READ_SUPPORT"])
    results_df.to_csv(output_file, index=False)
    print(f"Phased SV assignments written to: {output_file}")

haplotypes = get_read_haplotypes(bam_file)
print(f"Haplotypes extracted for {len(haplotypes)} chromosomes ✅")

assign_sv_to_haplotype(sv_vcf, haplotypes, output_file, min_support=5)