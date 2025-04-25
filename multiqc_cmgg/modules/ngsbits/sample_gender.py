#!/usr/bin/env python
import pysam
from cyvcf2 import VCF

#distribution on x and y chromosome (chat gpt)
def gender_xy(bam_file,max_female=0.06,min_male=0.09, include_single_end_reads=False):
    samfile = pysam.AlignmentFile(bam_file, "rb")

    chr_x_reads = 0
    chr_y_reads = 0

    for read in samfile.fetch():
        if read.is_unmapped or read.is_duplicate or read.is_secondary or read.is_supplementary:
            continue

        if not include_single_end_reads and (not read.is_paired or not read.is_proper_pair):
            continue

        ref_name = samfile.get_reference_name(read.reference_id)

        if ref_name in ["X", "chrX"]:
            chr_x_reads += 1
        elif ref_name in ["Y", "chrY"]:
            chr_y_reads += 1

    samfile.close()

    y_x_ratio = chr_y_reads / chr_x_reads if chr_x_reads > 0 else 0.0

    if y_x_ratio <= max_female:
        gender = "female"
    elif y_x_ratio >= min_male:
        gender = "male"
    else:
        gender = "unknown"

    return {
        "gender": gender,
        "add_info": {
            "reads_chrX": chr_x_reads,
            "reads_chrY": chr_y_reads,
            "y_x_ratio": round(y_x_ratio, 4),
        }
    }

#fraction of heterozygous variants, Use of build file ?, BAM to vcf with bcftools
def gender_hetx(vcf_file, max_male=0.05, min_female=0.25):
    het = 0
    total = 0

    for variant in VCF(vcf_file):
        if variant.CHROM not in ["X", "chrX"]:
            continue
        if variant.num_alt == 0:
            continue  # Skip reference-only sites

        # For diploid samples only
        gt = variant.genotypes[0]  # e.g. [0, 1, True]
        if gt[0] is None or gt[1] is None:
            continue

        total += 1
        if gt[0] != gt[1]:
            het += 1

    het_fraction = het / total if total > 0 else 0.0

    if het_fraction >= min_female:
        gender = "female"
    elif het_fraction <= max_male:
        gender = "male"
    else:
        gender = "unknown"

    return {
        "gender": gender,
        "add_info": {
            "total_snps": total,
            "heterozygous_snps": het,
            "het_fraction": round(het_fraction, 4),
        }
    }


#coverage of SRY gene

def gender_sry(bam_file, build="hg38", min_coverage=20.0):
    # Define SRY region depending on build (hardcoded, softcoding?)
    if build == "hg19":
        chrom = "Y"
        start = 2781507
        end = 2781984
    elif build == "hg38":
        chrom = "Y"
        start = 2782183
        end = 2782661
    else:
        raise ValueError("Unsupported genome build. Use 'hg19' or 'hg38'.")

    samfile = pysam.AlignmentFile(bam_file, "rb")
    
    coverage_total = 0
    base_count = end - start + 1

    try:
        for pileupcolumn in samfile.pileup(chrom, start, end, truncate=True):
            if start <= pileupcolumn.pos <= end:
                coverage_total += pileupcolumn.nsegments
    except ValueError:
        # Region not found in BAM (e.g., no chrY present at all)
        base_count = 1  # Prevent divide-by-zero
        coverage_total = 0

    samfile.close()

    avg_coverage = coverage_total / base_count

    gender = "male" if avg_coverage >= min_coverage else "female"

    return {
        "gender": gender,
        "add_info": {
            "avg_sry_coverage": round(avg_coverage, 2),
            "sry_region": f"{chrom}:{start}-{end}",
            "min_male_cov": min_coverage
        }
    }