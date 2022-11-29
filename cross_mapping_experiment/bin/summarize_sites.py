import json
import os
from collections import defaultdict
from operator import itemgetter
import multiprocessing
from math import ceil
import argparse

from iggtools.common.utils import tsprint, command, InputStream, OutputStream, select_from_tsv


DEFAULT_SITE_TYPES = ["ref", "snp"]


snps_pileup_schema = {
    "ref_id": str,
    "ref_pos": int,
    "ref_allele": str,
    "depth": int,
    "count_a": int,
    "count_c": int,
    "count_g": int,
    "count_t": int,
    "major_allele": str,
    "minor_allele": str,
    "major_allele_freq": float,
    "minor_allele_freq": float,
    "allele_counts": int,
}


true_sites_schema = {
    "refid": str,
    "qryid": str,
    "pos1": str,
    "pos2": str,
    "refallele": str,
    "altallele": str,
    "sitetype": str,
}


sites_stats_schema = {
    "refid": str,
    "sitetype": str,
    "sites_total": int,
    "sites_uncalled": int,
    "sites_correct": int,
    "sites_incorrect": int,
}
# sites_unaligned = total - uncalled - correct - incorrect !!


def format_data(x, decimal = ".3f"):
    return format(x, decimal) if isinstance(x, float) else str(x)


def get_site_id(ref_id, ref_pos, ref_allele):
    return f"{ref_id}|{ref_pos}|{ref_allele}"


def gen_myid(ref_id, site_type):
    return f"{ref_id}|{site_type}"


def read_true_sites(truesites_fp):
    # Read-only global variables
    true_sites_dict = defaultdict(dict)
    sites_summary_dict = defaultdict(dict)
    with InputStream(truesites_fp) as stream:
        for row in select_from_tsv(stream, selected_columns=true_sites_schema, result_structure=dict):
            ref_id = row["refid"]
            ref_pos = str(int(float(row["pos1"])))
            site_id = get_site_id(ref_id, ref_pos, row["refallele"])
            true_sites_dict[site_id] = row
            site_type = row["sitetype"]
            my_id = gen_myid(ref_id, site_type)
            if my_id in sites_summary_dict:
                sites_summary_dict[my_id]["sites_total"] += 1
            else:
                sites_stats = {
                    "refid": ref_id,
                    "sitetype":site_type,
                    "sites_total": 1,
                    "sites_uncalled": 0,
                    "sites_correct": 0,
                    "sites_incorrect": 0,
                }
                sites_summary_dict[my_id] = sites_stats
    return (true_sites_dict, sites_summary_dict)


def process_pileup(pileup_fp, truesites_fp, out_fp):
    true_sites_dict, sites_summary_dict = read_true_sites(truesites_fp)

    with OutputStream(out_fp) as ostream:
        with InputStream(pileup_fp) as stream:
            for row in select_from_tsv(stream, selected_columns=snps_pileup_schema, result_structure=dict):
                ref_id = row["ref_id"]
                ref_pos = row["ref_pos"]
                site_id = get_site_id(ref_id, ref_pos, row["ref_allele"])
                if site_id in true_sites_dict:
                    site_type = true_sites_dict[site_id]["sitetype"]
                    my_id =gen_myid(ref_id, site_type)

                    if row["major_allele_freq"] == row["minor_allele_freq"]:
                        sites_summary_dict[my_id]["sites_uncalled"] += 1
                    else:
                        alt_allele = true_sites_dict[site_id]["altallele"]
                        major_allele = row["major_allele"]
                        major_allele_readcounts = round(row["depth"] * row["major_allele_freq"])

                        if major_allele == alt_allele:
                            sites_summary_dict[my_id]["sites_correct"] += 1
                            curr_correct = "c"
                        else:
                            sites_summary_dict[my_id]["sites_incorrect"] += 1
                            curr_correct = "i"
                        ostream.write("\t".join(map(format_data, [ref_id, ref_pos, row["ref_allele"], major_allele, major_allele_readcounts, curr_correct])) + "\n")
    return sites_summary_dict


def main():
    p = argparse.ArgumentParser(prog="python summarize_pileup.py", description='compute sites level summary per pileup file')
    p.add_argument(
        "--pileup_fp", required=True,
        type=str,
        help=f"Path of pileup file.")
    p.add_argument(
        "--truesites_fp", required=True,
        type=str,
        help=f"Path to reference set files.")
    p.add_argument(
        "--summary_fp", required=True,
        type=str,
        help=f"Path to output summary file.")
    p.add_argument(
        "--out_fp", required=True,
        type=str,
        help=f"Path to output major allele file.")

    args = p.parse_args()
    if os.path.exists(args.pileup_fp):
        sites_summary_dict = process_pileup(args.pileup_fp, args.truesites_fp, args.out_fp)
    with OutputStream(args.summary_fp) as stream:
        stream.write("\t".join(sites_stats_schema.keys()) + "\n")
        if os.path.exists(args.pileup_fp):
            for _, record in sites_summary_dict.items() :
                stream.write("\t".join(map(str, record.values())) + "\n")

main()
