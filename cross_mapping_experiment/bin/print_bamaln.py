
import argparse
import os
import sys
import configparser
from collections import defaultdict
from Bio import SeqIO
from math import ceil
import pysam
import numpy as np


def format_data(x, decimal = ".3f"):
  return format(x, decimal) if isinstance(x, float) else str(x)


def print_aln_record(aln):
  """ Print one alignment record for one single read """
  row = [aln.reference_name, aln.reference_start, aln.reference_end, aln.query_name,
        aln.query_alignment_start, aln.query_alignment_end,
        "R1" if aln.is_read1 else "R2", aln.is_secondary, aln.is_reverse, aln.mate_is_reverse,
        aln.is_paired, aln.is_proper_pair, aln.next_reference_name,
        aln.template_length]
  align_len = len(aln.query_alignment_sequence)
  query_len = aln.query_length
  readq = "{:.2f}".format(np.mean(aln.query_qualities))
  mapq = aln.mapping_quality
  mismatches = dict(aln.tags)['NM']
  mapid = "{:.2f}".format(100 * (align_len - mismatches) / float(align_len))
  alncov = "{:.2f}".format(align_len / float(query_len))
  row = row + [align_len, query_len, readq, mapq, mismatches, mapid, alncov]
  return row



def main():
    p = argparse.ArgumentParser(prog="python print_bamaln.py", description='compute sites level summary per pileup file')
    p.add_argument(
        "--species_id", required=True,
        type=str,
        help=f"Path to reference set files.")
    p.add_argument(
        "--on_bin", required=True,
        type=str,
        help=f"Path to reference set files.")
    p.add_argument(
        "--off_bin", required=True,
        type=str,
        help=f"Path to reference set files.")

    args = p.parse_args()
    species_id = args.species_id
    on_bin = args.on_bin
    off_bin = args.off_bin

    basedir = "/mnt/cz/20220602_close_species/"
    bamfile = f"{basedir}/{species_id}/{on_bin}/6_midas_exp1/{off_bin}/cov_20X/temp/snps/repgenomes.bam"
    outfile = f"{basedir}/{species_id}/bamalns/{on_bin}/{off_bin}.tsv"
    outdir = f"{basedir}/{species_id}/bamalns/{on_bin}"

    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    with open(outfile, "w") as ofile:
        with pysam.AlignmentFile(bamfile) as infile:
            for aln in infile:
                row = print_aln_record(aln)
                ofile.write("\t".join(map(format_data, row)) + "\n")

main()
