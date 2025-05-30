#!/usr/bin/env python3

import argparse


def read_bed(bed_file):
    """
    Read a BED file and store its contents in a dictionary.

    Parameters:
    - bed_file (str): Path to the BED file.

    Returns:
    - dict: A dictionary where keys are chromosome names and values are lists of tuples.
            Each tuple contains (start, end, gene, strand).
    """
    bed_dict = {}
    with open(bed_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 6:
                print(f"Warning: Line '{line.strip()}' in BED file has fewer than 6 columns, skipping.")
                continue
            chrom, start, end, gene, score, strand = fields[:6]
            if chrom not in bed_dict:
                bed_dict[chrom] = []
            bed_dict[chrom].append((int(start), int(end), gene, strand))
    return bed_dict


def process_vcf(vcf_file, bed_dict, output_file):
    """
    Process a VCF file by adding gene and strand information from a BED dictionary,
    and write the result to an output file.

    Parameters:
    - vcf_file (str): Path to the input VCF file.
    - bed_dict (dict): Dictionary containing BED data.
    - output_file (str): Path to the output VCF file with added gene and strand info.
    """
    with open(vcf_file, 'r') as vcf, open(output_file, 'w') as out:
        # Copy the header line directly and add GENE and STRAND column names
        header = next(vcf).strip()
        out.write(header + '\tGENE\tSTRAND\n')

        for line in vcf:
            fields = line.strip().split('\t')
            chrom, pos = fields[0], int(fields[1])
            gene_name = ''
            strand = ''
            if chrom in bed_dict:
                for start, end, gene, s in bed_dict[chrom]:
                    if start <= pos <= end:
                        gene_name = gene
                        strand = s
                        break
            # Append gene name and strand information to the end of the line
            out.write('\t'.join(fields) + f'\t{gene_name}\t{strand}\n')


def main():
    parser = argparse.ArgumentParser(
        description='This script adds gene and strand information from a BED file to a VCF-like file.',
        epilog='Example usage: python add_gene_strand_to_vcf.py -i input.vcf -r regions.bed -o output.vcf'
    )
    parser.add_argument('-i', '--input_vcf', required=True, help='Path to the input VCF-like file.')
    parser.add_argument('-r', '--input_bed', required=True, help='Path to the input BED file.')
    parser.add_argument('-o', '--output_vcf', required=True, help='Path to the output VCF-like file.')

    args = parser.parse_args()

    # Read BED file
    bed_dict = read_bed(args.input_bed)
    # Process VCF file
    process_vcf(args.input_vcf, bed_dict, args.output_vcf)


if __name__ == '__main__':
    main()