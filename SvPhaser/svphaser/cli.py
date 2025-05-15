import argparse
import os
from multiprocessing import Pool
import pandas as pd
import glob
import pysam
from collections import defaultdict, Counter


def get_read_haplotypes(bam_file, chrom):
    haplotypes = defaultdict(list)
    bam = pysam.AlignmentFile(bam_file, "rb")
    for read in bam.fetch(chrom):
        if read.has_tag("HP"):
            hp_tag = read.get_tag("HP")
            haplotypes[read.reference_name].append((read.reference_start, read.reference_end, hp_tag))
    bam.close()
    return haplotypes


def assign_sv_to_haplotype_by_chromosome(args):
    chrom, bam_file, sv_vcf, output_folder, min_support = args
    haplotypes = get_read_haplotypes(bam_file, chrom)

    compression = 'gzip' if sv_vcf.endswith('.gz') else None

    sv_df = pd.read_csv(
        sv_vcf,
        comment="#",
        sep="\t",
        header=None,
        compression=compression
    )

    sv_df.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]
    sv_df = sv_df[sv_df["CHROM"] == chrom]

    results = []

    for _, row in sv_df.iterrows():
        pos = int(row["POS"])
        sv_start = pos
        sv_end = pos + 1

        haps = [
            hp for start, end, hp in haplotypes.get(chrom, [])
            if start <= sv_start <= end
        ]

        support_count = len(haps)
        assigned_hap = "unphased"

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

    os.makedirs(output_folder, exist_ok=True)
    results_df = pd.DataFrame(results, columns=["CHROM", "POS", "SV_ID", "HAPLOTYPE", "READ_SUPPORT"])
    results_df.to_csv(f"{output_folder}/SV_phasing_{chrom}.csv", index=False)
    print(f"{chrom} done ✅")


def merge_csvs(folder, output_file):
    csv_files = glob.glob(f"{folder}/SV_phasing_chr*.csv")
    df_list = [pd.read_csv(file) for file in csv_files]
    merged_df = pd.concat(df_list, ignore_index=True)
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    merged_df.to_csv(output_file, index=False)


def write_phased_vcf(csv_file, vcf_input, vcf_output):
    df = pd.read_csv(csv_file)
    phasing_dict = {
        (row["CHROM"], row["POS"]): (row["HAPLOTYPE"], row["READ_SUPPORT"])
        for _, row in df.iterrows()
    }

    with open(vcf_input, "r") as fin, open(vcf_output, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue

            parts = line.strip().split("\t")
            chrom = parts[0]
            pos = int(parts[1])

            key = (chrom, pos)
            haplotype, read_support = phasing_dict.get(key, ("unphased", 0))

            if haplotype == "HP1":
                gt = "1|0"
            elif haplotype == "HP2":
                gt = "0|1"
            elif haplotype == "HP1_HP2":
                gt = "1|1"
            else:
                gt = "./."

            parts[8] = "GT:DP"
            parts[9] = f"{gt}:{read_support}"

            fout.write("\t".join(parts) + "\n")


def main():
    parser = argparse.ArgumentParser(description="SvPhaser: Phases SVs from a phased BAM + unphased SV VCF.")
    parser.add_argument("--phased_bam", required=True, help="Path to phased BAM file")
    parser.add_argument("--unphased_vcf", required=True, help="Path to unphased SV VCF file")
    parser.add_argument("--output", required=True, help="Output folder")
    parser.add_argument("--min_support", type=int, default=10, help="Minimum read support (default=10)")

    args = parser.parse_args()

    chrom_folder = os.path.join(args.output, "chromosome_csvs")
    merged_folder = os.path.join(args.output, "merged")

    chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    tasks = [(chrom, args.phased_bam, args.unphased_vcf, chrom_folder, args.min_support) for chrom in chromosomes]

    with Pool() as pool:
        pool.map(assign_sv_to_haplotype_by_chromosome, tasks)

    merged_csv = os.path.join(merged_folder, "SV_phasing_full.csv")
    merge_csvs(chrom_folder, merged_csv)

    input_name = os.path.splitext(os.path.basename(args.unphased_vcf))[0]
    phased_vcf = os.path.join(merged_folder, f"{input_name}_phased.vcf")
    write_phased_vcf(merged_csv, args.unphased_vcf, phased_vcf)

    # Delete temporary chromosome CSVs
    import glob
    csv_files = glob.glob(os.path.join(chrom_folder, "SV_phasing_chr*.csv"))
    for file in csv_files:
        os.remove(file)
    print("Temporary chromosome CSVs deleted ✅")

    print("SvPhaser complete ✅")


if __name__ == "__main__":
    main()
