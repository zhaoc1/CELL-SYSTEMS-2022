import pysam
import numpy as np
import argparse
from collections import defaultdict

fastmode = True


def bt2_mapq_end2end(AS, XS=None, scMin=-90.6):
    # http://biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html
    '''scMin = minScore'''
    if XS == None:
        XS = scMin-1
    if XS > AS:
        return None
    diff = abs(scMin)
    bestOver = AS-scMin
    bestdiff = abs(abs(AS)-abs(XS))
    if XS < scMin:
        if bestOver >= diff*0.8:
            return 42
        elif bestOver >= diff*0.7:
            return 40
        elif bestOver >= diff*0.61:
            return 24
        elif bestOver >= diff*0.5:
            return 23
        elif bestOver >= diff*0.42:
            return 8
        elif bestOver >= diff*0.3:
            return 3
        else:
            return 0
    else:
        if bestdiff >= diff*0.9:
            if bestOver == diff:
                return 39
            else:
                return 33
        elif bestdiff >= diff*0.8:
            if bestOver == diff:
                return 38
            else:
                return 27
        elif bestdiff >= diff*0.97:
            if bestOver == diff:
                return 37
            else:
                return 26
        elif bestdiff >= diff*0.6:
            if bestOver == diff:
                return 36
            else:
                return 22
        elif bestdiff >= diff*0.5:
            if bestOver == diff:
                return 35
            elif bestOver >= diff*0.84:
                return 25
            elif bestOver >= diff*0.68:
                return 16
            elif bestOver >= diff*0.68:
                return 5
        elif bestdiff >= diff*0.4:
            if bestOver == diff:
                return 34
            elif bestOver >= diff*0.84:
                return 21
            elif bestOver >= diff*0.68:
                return 14
            else:
                return 4
        elif bestdiff >= diff*0.3:
            if bestOver == diff:
                return 32
            elif bestOver >= diff*0.88:
                return 18
            elif bestOver >= diff*0.67:
                return 15
            else:
                return 3
        elif bestdiff >= diff*0.2:
            if bestOver == diff:
                return 31
            elif bestOver >= diff*0.88:
                return 17
            elif bestOver >= diff*0.67:
                return 11
            else:
                return 0
        elif bestdiff >= diff*0.1:
            if bestOver == diff:
                return 30
            elif bestOver >= diff*0.88:
                return 12
            elif bestOver >= diff*0.67:
                return 7
            else:
                return 0
        elif bestdiff > 0:
            if bestOver >= diff*0.67:
                return 6
            else:
                return 2
        else:
            if bestOver >= diff*0.68:
                return 1
            else:
                return 0


def hamming_distance(str1, str2):
    """ Compute the Hamming distance between two strings """
    assert len(str1) == len(str2), f"Two input strings for hamming_distance are different length."
    hd = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            hd += 1
    return hd


def compute_mapq(aln):
    """ for record """
    # This is the original single end reads case
    AS = dict(aln.tags)['AS']
    XS = dict(aln.tags)['XS'] if "XS" in dict(aln.tags) else None
    min_score = -0.6 + -0.6 * aln.query_length
    cal_mapq = bt2_mapq_end2end(AS, XS, min_score)


def check_proper_pair(aln):
    """ This is not continue anymore, just for record """
    if not aln.is_paired:
        return False
    # 2021-05-20: try to reproduce the is_proper_pair
    if not (aln.is_reverse ^ aln.mate_is_reverse):
        return False # reads aligned to different strand
    if aln.reference_name != aln.next_reference_name:
        return False # breads aligned to the same contig

    is_overlap = intervals_overlap((alns["fwd"].reference_start, alns["fwd"].reference_end), (alns["rev"].reference_start, alns["rev"].reference_end))
    # Can also just compute the overlap here
    if alns["rev"].reference_start >= alns["fwd"].reference_end:
        print(f"no overlap {is_overlap}")
    elif alns["rev"].reference_start >= alns["fwd"].reference_start:
        print(f"proper pair {is_overlap}")
    else:
        print(f"proper pair {is_overlap}")
    return True


def compute_single_read_mismatch(aln):
    # 2021-05-19: assert the compute of NM - DONE
    ref_seq = aln.get_reference_sequence()
    qry_seq = aln.query_alignment_sequence
    aligned_pos = aln.get_aligned_pairs()
    nm = 0
    ro = []
    qo = []
    for i in range(0,len(aligned_pos)):
        if aligned_pos[i][0] is None:
            qo.append("-")
        else:
            qo.append(qry_seq[aligned_pos[i][0] - aln.query_alignment_start])
            nm += 1

        if aligned_pos[i][1] is None:
            nm += 1
            ro.append("-")
        else:
            ro.append(ref_seq[aligned_pos[i][1] - aln.reference_start])

    ro = "".join(ro).upper()
    qo = "".join(qo).upper()
    hd_w_gaps = hamming_distance(ro, qo)
    assert (hd_w_gaps == dict(aln.tags)['NM']), "".join(ro) + "\n" + "".join(qo).upper()

    if False:
        if aligned_pos[i][0] is None or aligned_pos[i][1] is None:
            nm += 1
        else:
            ro.append(ref_seq[aligned_pos[i][1] - aln.reference_start])
            qo.append(qry_seq[aligned_pos[i][0] - aln.query_alignment_start])
        hd = hamming_distance("".join(ro).upper(), "".join(qo).upper())
        assert(nm+hd == dict(aln.tags)['NM']), "".join(ro) + "\n" + "".join(qo).upper()


def compute_single_nm(aln):
    # 2021-05-19: assert the compute of NM - DONE
    ref_seq = aln.get_reference_sequence()
    qry_seq = aln.query_alignment_sequence
    aligned_pos = aln.get_aligned_pairs()
    nm = 0
    ro = []
    qo = []
    for i in range(0,len(aligned_pos)):
        if aligned_pos[i][0] is None or aligned_pos[i][1] is None:
            nm += 1
        else:
            ro.append(ref_seq[aligned_pos[i][1] - aln.reference_start])
            qo.append(qry_seq[aligned_pos[i][0] - aln.query_alignment_start])
        hd = hamming_distance("".join(ro).upper(), "".join(qo).upper())
        assert(nm+hd == dict(aln.tags)['NM']), "".join(ro) + "\n" + "".join(qo).upper()
    return (nm+hd)


def interval(a, b):
    return (min(a, b), max(a, b))


def intervals_overlap(p, q):
    return max(0.0,  min(p[1], q[1]) - max(p[0], q[0]) + 1)


def compute_mismatches_before_overlap(aln, reads_overlap):
    ref_seq = aln.get_reference_sequence() # reference sequence that is covered by the alignment of the read to the reference.
    qry_seq = aln.query_alignment_sequence # aligned portion of the read.
    aligned_pos = aln.get_aligned_pairs()

    nm = 0 #mismatches
    ro = []
    qo = []
    for i in range(0,len(aligned_pos)):
        if aligned_pos[i][0] is not None and aligned_pos[i][0] >= aln.query_alignment_end - reads_overlap:
            break

        if aligned_pos[i][0] is None:
            qo.append("-")
            nm += 1
        else:
            qo.append(qry_seq[aligned_pos[i][0] - aln.query_alignment_start])

        if aligned_pos[i][1] is None:
            ro.append("-")
            nm += 1
        else:
            ro.append(ref_seq[aligned_pos[i][1] - aln.reference_start])

    ro = "".join(ro).upper()
    qo = "".join(qo).upper()
    hd_w_gaps = hamming_distance(ro, qo)

    alnlen_wo_gaps = len([_ for _ in qo if _ is not "-"])
    assert alnlen_wo_gaps == aln.query_alignment_length - reads_overlap, f"wrong alignment length"


    a = []
    b = []
    nm = 0
    for i in range(0,len(aligned_pos)):
        if aligned_pos[i][0] is not None and aligned_pos[i][0] >= aln.query_alignment_end - reads_overlap:
            break

        if aligned_pos[i][0] is None or aligned_pos[i][1] is None:
            nm += 1
        else:
            b.append(qry_seq[aligned_pos[i][0] - aln.query_alignment_start])
            a.append(ref_seq[aligned_pos[i][1] - aln.reference_start])
    hd_wo_gaps = hamming_distance("".join(a).upper(), "".join(b).upper())

    assert hd_w_gaps == hd_wo_gaps + nm
    ## As of now, I can confirm the gaps are treated as aligned mismatches ... This might not be expected for alncov. But lets worry about this later.

    return hd_w_gaps


def compute_mismatches_forward(aln, reads_overlap):
    ref_seq = aln.get_reference_sequence() # reference sequence that is covered by reads alignment
    qry_seq = aln.query_alignment_sequence # aligned portion of the read

    aligned_pos = aln.get_aligned_pairs()

    ngaps_o = 0 #mismatches
    ngaps_i = 0 #mismatches
    ro = []
    qo = []
    ri = []
    qi = []

    for i in range(0,len(aligned_pos)):
        if aligned_pos[i][0] is not None and aligned_pos[i][0] >= aln.query_alignment_end - reads_overlap:
            if aligned_pos[i][0] is None:
                qi.append("-")
                ngaps_i += 1
            else:
                qi.append(qry_seq[aligned_pos[i][0] - aln.query_alignment_start])

            if aligned_pos[i][1] is None:
                ri.append("-")
                ngaps_i += 1
            else:
                ri.append(ref_seq[aligned_pos[i][1] - aln.reference_start])
        else:
            # Outside the overlap region
            if aligned_pos[i][0] is None:
                qo.append("-")
                ngaps_o += 1
            else:
                qo.append(qry_seq[aligned_pos[i][0] - aln.query_alignment_start])

            if aligned_pos[i][1] is None:
                ro.append("-")
                ngaps_o += 1
            else:
                ro.append(ref_seq[aligned_pos[i][1] - aln.reference_start])

    ro = "".join(ro)
    qo = "".join(qo)
    hd_out_w_gaps = hamming_distance(ro.upper(), qo.upper())

    ri = "".join(ri)
    qi = "".join(qi)
    hd_in_w_gaps = hamming_distance(ri.upper(), qi.upper())

    return (ro, qo, hd_out_w_gaps, ri, qi, hd_in_w_gaps)


def compute_mismatches_reverse(aln, reads_overlap):
    ref_seq = aln.get_reference_sequence() # reference sequence that is covered by reads alignment
    qry_seq = aln.query_alignment_sequence # aligned portion of the read

    aligned_pos = aln.get_aligned_pairs()

    ngaps_o = 0 #mismatches
    ngaps_i = 0 #mismatches
    ro = []
    qo = []
    ri = []
    qi = []

    for i in range(0,len(aligned_pos)):
        if aligned_pos[i][0] is not None and aligned_pos[i][0] <= aln.query_alignment_start + reads_overlap - 1:
            if aligned_pos[i][0] is None:
                qi.append("-")
                ngaps_i += 1
            else:
                qi.append(qry_seq[aligned_pos[i][0] - aln.query_alignment_start])

            if aligned_pos[i][1] is None:
                ri.append("-")
                ngaps_i += 1
            else:
                ri.append(ref_seq[aligned_pos[i][1] - aln.reference_start])
        else:
            # Outside the overlap region
            if aligned_pos[i][0] is None:
                qo.append("-")
                ngaps_o += 1
            else:
                qo.append(qry_seq[aligned_pos[i][0] - aln.query_alignment_start])

            if aligned_pos[i][1] is None:
                ro.append("-")
                ngaps_o += 1
            else:
                ro.append(ref_seq[aligned_pos[i][1] - aln.reference_start])

    ro = "".join(ro)
    qo = "".join(qo)
    hd_out_w_gaps = hamming_distance(ro.upper(), qo.upper())

    ri = "".join(ri)
    qi = "".join(qi)
    hd_in_w_gaps = hamming_distance(ri.upper(), qi.upper())

    return (ro, qo, hd_out_w_gaps, ri, qi, hd_in_w_gaps)


def main_20210521_prototype():
    # Complete on 2021-05-21. We keep the REV reads and correct the FWD reads.
    p = argparse.ArgumentParser(prog="python summarize_pileup.py", description='compute sites level summary per pileup file')
    p.add_argument(
        "--bamfile", required=True,
        type=str,
        help=f"Path of pileup file.")
    p.add_argument(
        "--outfile", required=True,
        type=str,
        help=f"Path to reference set files.")

    args = p.parse_args()

    outfile = args.outfile
    bamfile = args.bamfile


    alndict = defaultdict(dict)
    tally = 0
    #defaultdict(lambda: defaultdict(int))

    with pysam.AlignmentFile(bamfile) as infile:
        for aln in infile:
            if aln.is_secondary:
                continue
            if not aln.is_proper_pair:
                continue
            if aln.is_reverse:
                alndict[aln.query_name]["rev"] = aln
            else:
                alndict[aln.query_name]["fwd"] = aln

            tally +=1
            if tally == 20000:
                break

    for query_name, alns in alndict.items():
        if len(alns) != 2:
            continue

        if dict(alns["rev"].tags)['NM'] < 4:
            continue

        # Common attributes
        readq = "{:.3f}".format(np.mean(alns["fwd"].query_qualities + alns["rev"].query_qualities))
        mapq = max(alns["fwd"].mapping_quality, alns["rev"].mapping_quality)
        frag_len = abs(alns["fwd"].template_length) # number of bases from the left most mapped base to the rightmost mapped base on the reference.

        # The alignment coverage should not be effected by overlap.
        # However, we should double check whethe counts gaps as aligned ?!
        align_len = alns["fwd"].query_alignment_length + alns["rev"].query_alignment_length
        query_len = alns["fwd"].query_length + alns["rev"].query_length
        alncov = "{:.3f}".format(align_len / float(query_len))

        # Forwad strand is not necessarily the R1 reads ...
        reads_overlap = intervals_overlap((alns["fwd"].reference_start, alns["fwd"].reference_end - 1), (alns["rev"].reference_start, alns["rev"].reference_end - 1))

        if reads_overlap:
            # Record REV strand
            # For the FWD read, we read in the aligned sequences until overlap
            aln = alns["fwd"]

            #hd_w_gaps = compute_mismatches_before_overlap(aln, reads_overlap)
            (ro, qo, hd_out_w_gaps, ri, qi, hd_in_w_gaps) = compute_mismatches_forward(aln, reads_overlap)
            assert hd_out_w_gaps + hd_in_w_gaps == dict(alns["fwd"].tags)['NM']

            mismatches = dict(alns["rev"].tags)['NM'] + hd_out_w_gaps

            #length of the aligned query sequence. #<==== TODO: I think this one should not include the deletion
            align_len = alns["rev"].query_alignment_length + aln.query_alignment_length - reads_overlap
            mapid = "{:.3f}".format(100 * (align_len - mismatches) / float(align_len))

            alnf = alns["fwd"]
            alnr = alns["rev"]


            row = [f"{alnf.reference_name}:{alnf.reference_start}-{alnf.reference_end}, refCoveredAln: {alnf.get_overlap(aln.reference_start, alnf.reference_end)}",
                    "R1" if alnf.is_read1 else "R2",
                    alnf.get_reference_sequence(),
                    alnf.get_aligned_pairs(),
                    f"{alnf.query_name}:{alnf.query_alignment_start}-{alnf.query_alignment_end}, alnLen:{alnf.query_alignment_length}, query length:{alnf.query_length}",
                    alnf.query_alignment_sequence,
                    f"overlap between reads:{reads_overlap},  mismatches for FWD before overlap:{hd_out_w_gaps}, mismatches withint overlap:{hd_in_w_gaps}, total mismatches:{mismatches}, {dict(alnf.tags)['NM']}, {dict(alnr.tags)['NM']}",
                    f"{ro}\n{qo}\n{ri}\n{qi}",
                    f"aligned_length:{align_len}, template length:{frag_len}",
                    f"{alnr.reference_name}:{alnr.reference_start}-{alnr.reference_end}, refCoveredAln: {alnr.get_overlap(alnr.reference_start, alnr.reference_end)}",
                    "R1" if alnr.is_read1 else "R2",
                    alnr.get_reference_sequence(),
                    alnr.get_aligned_pairs(),
                    f"{alnr.query_name}:{alnf.query_alignment_start}-{alnr.query_alignment_end}, alnLen:{alnr.query_alignment_length}, query length:{alnr.query_length}",
                    alnr.query_alignment_sequence,
                    ]


            print("\n".join(map(str, row)) + "\n-------------------------------------------------------------------------------------------------\n")
            if False:
                if alnf.query_alignment_sequence[-reads_overlap:] != alnr.query_alignment_sequence[:reads_overlap]:

                        print(alnf.query_alignment_sequence)
                        print(alnr.query_alignment_sequence)
                        print(alnf.query_alignment_sequence[-reads_overlap:])
                        print(alnr.query_alignment_sequence[:reads_overlap])
                        assert alnf.query_alignment_sequence[-reads_overlap:] == alnr.query_alignment_sequence[:reads_overlap], f"the overlap region has mismatch"
        else:
            mismatches = dict(alns["fwd"].tags)['NM'] + dict(alns["rev"].tags)['NM']
            mapid = "{:.3f}".format(100 * (align_len - mismatches) / float(align_len))


    if False:
        # If there are no overlaps between two reads, then we can simply add up
        if aln.query_name not in alndict:
            if aln.is_read1:
                alndict[aln.query_name] = {"mapq": mapq, "edist": edist, "alnlen": alnlen, "qrylen": qrylen}
            else:
                # if the same query
                mapq = max(mapq, aln.mapping_quality)
                edist += edist + dict(aln.tags)['NM']
                alnlen += alnlen + aln.query_alignment_length ## This is not necessarily true. Because there would be overlap reads ...
                qrylen += aln.query_length
            #in order to compute the mapid, you need to number: edit distance and alignment_length
            mapid = "{:.3f}".format(100 * (alnlen - edist) / alnlen)
            alncov = "{:.3f}".format(alnlen / qrylen)


def print_read_pairs(bamfile, outfile):

    outfile = args.outfile
    bamfile = args.bamfile


    alndict = defaultdict(dict)
    tally = 0

    with pysam.AlignmentFile(bamfile) as infile:
        for aln in infile:
            if aln.is_secondary:
                continue
            if not aln.is_proper_pair:
                continue
            if aln.is_reverse:
                alndict[aln.query_name]["rev"] = aln
            else:
                alndict[aln.query_name]["fwd"] = aln

            tally +=1
            if fastmode:
                if tally == 100000:
                    break

    with open(outfile, "w") as ofile:
        for query_name, alns in alndict.items():
            # Skip if not paired
            if len(alns) != 2:
                continue

            #if dict(alns["rev"].tags)['NM'] < 4:
            #    continue

            # Common attributes
            readq = "{:.3f}".format(np.mean(alns["fwd"].query_qualities + alns["rev"].query_qualities))
            mapq = max(alns["fwd"].mapping_quality, alns["rev"].mapping_quality)
            frag_len = abs(alns["fwd"].template_length) # number of bases from the left most mapped base to the rightmost mapped base on the reference.

            # The alignment coverage should not be effected by overlap.
            # However, we should double check whethe counts gaps as aligned ?!
            align_len = alns["fwd"].query_alignment_length + alns["rev"].query_alignment_length
            query_len = alns["fwd"].query_length + alns["rev"].query_length
            alncov = "{:.3f}".format(align_len / float(query_len))

            # Forwad strand is not necessarily the R1 reads ...
            reads_overlap = intervals_overlap((alns["fwd"].reference_start, alns["fwd"].reference_end - 1), (alns["rev"].reference_start, alns["rev"].reference_end - 1))

            if reads_overlap:
                # Keep the REV read, split the FWD reads
                (ro, qo, hd_out_w_gaps, ri, qi, hd_in_w_gaps) = compute_mismatches_forward(alns["fwd"], reads_overlap)
                #assert hd_out_w_gaps + hd_in_w_gaps == dict(alns["fwd"].tags)['NM']

                # Keep the FWD read, split the REV reads
                (ro2, qo2, hd_out_w_gaps2, ri2, qi2, hd_in_w_gaps2) = compute_mismatches_reverse(alns["rev"], reads_overlap)
                #assert hd_out_w_gaps2 + hd_in_w_gaps2 == dict(alns["rev"].tags)['NM']

                # We only consider proper aligned pairsself
                # On top of that,
                if abs(hd_in_w_gaps - hd_in_w_gaps2) > 1:
                    continue

                mismatches = dict(alns["rev"].tags)['NM'] + hd_out_w_gaps
                mismatches2 = dict(alns["fwd"].tags)['NM'] + hd_out_w_gaps2

                #length of the aligned query sequence. #<==== TODO: I think this one should not include the deletion
                align_len = alns["rev"].query_alignment_length + aln.query_alignment_length - reads_overlap
                mapid = "{:.3f}".format(100 * (align_len - mismatches) / float(align_len))

                alnf = alns["fwd"]
                alnr = alns["rev"]

                row = [alnf.query_name, alnf.reference_name, alnf.reference_start, alnf.reference_end,
                        "FWD-R1" if alnf.is_read1 else "FWD-R2", dict(alnf.tags)['NM'],
                        alnf.query_length, alnr.query_length,
                        "REV-R1" if alnr.is_read1 else "REV-R2", dict(alnr.tags)['NM'],
                        reads_overlap,
                        hd_out_w_gaps, hd_in_w_gaps, mismatches,
                        hd_out_w_gaps2, hd_in_w_gaps2, mismatches2,
                        frag_len, align_len]

                ofile.write("\t".join(map(str, row)) + "\n")

                if False:
                    # For the purpose of debug whether the overlap fragment between read pairs are the same
                    #assert hd_in_w_gaps == hd_in_w_gaps2, f"Let's come back here after the QC"
                    if alnf.query_alignment_sequence[-reads_overlap:] != alnr.query_alignment_sequence[:reads_overlap]:
                            print(alnf.query_alignment_sequence)
                            print(alnr.query_alignment_sequence)
                            print(alnf.query_alignment_sequence[-reads_overlap:])
                            print(alnr.query_alignment_sequence[:reads_overlap])
                            assert alnf.query_alignment_sequence[-reads_overlap:] == alnr.query_alignment_sequence[:reads_overlap]

                if False:
                    # Print for Debug purpose
                    row = [f"{alnf.reference_name}:{alnf.reference_start}-{alnf.reference_end}, refCoveredAln: {alnf.get_overlap(aln.reference_start, alnf.reference_end)}",
                            "R1" if alnf.is_read1 else "R2",
                            alnf.get_reference_sequence(),
                            alnf.get_aligned_pairs(),
                            f"{alnf.query_name}:{alnf.query_alignment_start}-{alnf.query_alignment_end}, alnLen:{alnf.query_alignment_length}, query length:{alnf.query_length}",
                            alnf.query_alignment_sequence,
                            f"overlap between reads:{reads_overlap},  mismatches for FWD before overlap:{hd_out_w_gaps}, mismatches within overlap:{hd_in_w_gaps}, total mismatches:{mismatches}, {dict(alnf.tags)['NM']}, {dict(alnr.tags)['NM']}",
                            f"overlap between reads:{reads_overlap},  mismatches for REV after overlap:{hd_out_w_gaps2}, mismatches within overlap:{hd_in_w_gaps2}, total mismatches:{mismatches2}, {dict(alnf.tags)['NM']}, {dict(alnr.tags)['NM']}",
                            f"{ro}\n{qo}\n{ri}\n{qi}",
                            f"aligned_length:{align_len}, template length:{frag_len}",
                            f"{alnr.reference_name}:{alnr.reference_start}-{alnr.reference_end}, refCoveredAln: {alnr.get_overlap(alnr.reference_start, alnr.reference_end)}",
                            "R1" if alnr.is_read1 else "R2",
                            alnr.get_reference_sequence(),
                            alnr.get_aligned_pairs(),
                            f"{alnr.query_name}:{alnf.query_alignment_start}-{alnr.query_alignment_end}, alnLen:{alnr.query_alignment_length}, query length:{alnr.query_length}",
                            alnr.query_alignment_sequence,
                        ]
                    print("\n".join(map(str, row)) + "\n-------------------------------------------------------------------------------------------------\n")
            else:
                mismatches = dict(alns["fwd"].tags)['NM'] + dict(alns["rev"].tags)['NM']
                mapid = "{:.3f}".format(100 * (align_len - mismatches) / float(align_len))


def main():

    p = argparse.ArgumentParser(prog="python summarize_pileup.py", description='compute sites level summary per pileup file')
    p.add_argument(
        "--bamfile", required=True,
        type=str,
        help=f"Path of pileup file.")
    p.add_argument(
        "--outread", required=True,
        type=str,
        help=f"Path to reference set files.")
    p.add_argument(
        "--outpair", required=True,
        type=str,
        help=f"Path to reference set files.")

    args = p.parse_args()

    outfile1 = args.outread
    bamfile = args.bamfile
    outfile2 = args.outpair

    print_read(bamfile, outfile1)
    print_read_pairs(bamfile, outfile2)


def print_read(bamfile, outfile):
    tally = 0
    with open(outfile, "w") as ofile:
        with pysam.AlignmentFile(bamfile) as infile:
            alndict = defaultdict(dict)
            for aln in infile:
                if aln.is_paired:
                    row = [aln.reference_name, aln.reference_start, aln.reference_end,
                            aln.query_name, aln.query_alignment_start, aln.query_alignment_end,
                            aln.template_length, aln.mapping_quality, dict(aln.tags)['NM'],
                            "R1" if aln.is_read1 else "R2", "rev" if aln.is_reverse else "fwd",
                            aln.is_proper_pair, aln.is_reverse, aln.mate_is_reverse, aln.reference_name, aln.next_reference_name]

                    if False:
                        aligned_pos = aln.get_aligned_pairs()
                        alnlen = len(aligned_pos)
                        gapq = sum(aligned_pos[i][0] is None for i in range(0,len(aligned_pos)))
                        gapr = sum(aligned_pos[i][1] is None for i in range(0,len(aligned_pos)))
                        if gapq:
                            #assert aln.query_alignment_length == len(aligned_pos) - gapr - gapq, f"{aln.query_alignment_length}, {alnlen}, gapq"
                            tally += 1
                            print("\n".join(map(str, row)) + "\n-----------------------------------------------------\n")

                    #if not aln.is_proper_pair:
                    tally += 1
                    ofile.write("\t".join(map(str, row)) + "\n")
                    if fastmode:
                        if tally == 100000:
                            exit(0)


main()

#predefined tag NM means: Edit distance to the reference (number of changes
#necessary to make this equal the reference, excluding clipping)
