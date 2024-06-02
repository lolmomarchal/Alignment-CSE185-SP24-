import argparse
import time

from Bio import SeqIO

from ReferenceIndexing import referenceGenome as referenceGenome

# Parameters for seed and extend
DEFAULT_K = 10
DEFAULT_MISMATCH_THRESHOLD = 3
DEFAULT_GAP_PENALTY = -2
DEFAULT_MISMATCH_PENALTY = -1
DEFAULT_MATCH_SCORE = 2


def parse_arguments():
    parser = argparse.ArgumentParser(description="Seed and Extend BWT Aligner")
    parser.add_argument('-k', '--k', type=int, default=DEFAULT_K, help='Length of k-mer for seeding')
    parser.add_argument('-mt', '--mismatch_threshold', type=int, default=DEFAULT_MISMATCH_THRESHOLD,
                        help='Threshold for mismatches')
    parser.add_argument('-d', '--gap_penalty', type=int, default=DEFAULT_GAP_PENALTY, help='Penalty for gaps')
    parser.add_argument('-s', '--mismatch_penalty', type=int, default=DEFAULT_MISMATCH_PENALTY,
                        help='Penalty for mismatches')
    parser.add_argument('-m', '--match_score', type=int, default=DEFAULT_MATCH_SCORE, help='Score for matches')
    parser.add_argument('-r', '--reference_path', type=str, required=True, help='Path to reference genome file')
    parser.add_argument('-i', '--read_path', type=str, required=True, help='Path to reads file')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to output SAM file')
    return parser.parse_args()


def seed_and_extend_bwt(query, reference, k, mismatch_threshold):
    seeds = []
    for i in range(len(query) - k + 1):
        kmer = query[i:i + k]
        positions = reference.backward_search(kmer)
        for pos in positions:
            match_length, mismatches = extend_match(query, reference.sequence, i, pos, k, mismatch_threshold)
            if match_length > k:
                seeds.append((i, pos, match_length, mismatches))
    return seeds


def extend_match(query, reference, query_start, ref_start, seed_length, mismatch_threshold):
    match_length = seed_length
    mismatches = 0
    # Extend forward
    while (ref_start + match_length < len(reference) and
           query_start + match_length < len(query)):
        if reference[ref_start + match_length] != query[query_start + match_length]:
            mismatches += 1
            if mismatches > mismatch_threshold:
                break
        match_length += 1
    # Extend backward
    while (ref_start > 0 and query_start > 0):
        if reference[ref_start - 1] != query[query_start - 1]:
            mismatches += 1
            if mismatches > mismatch_threshold:
                break
        ref_start -= 1
        query_start -= 1
        match_length += 1
    return match_length, mismatches


def obtain_reads(path_reads):
    sequences = {}
    for record in SeqIO.parse(path_reads, "fastq"):
        sequences[record.id] = (str(record.seq), record.letter_annotations["phred_quality"])
    return sequences


def align(query, ref, query_start, ref_start, match_length, mismatch_penalty, gap_penalty, match_score,
          mismatch_threshold):
    alignment_score = 0
    mismatches = 0
    i, j = query_start, ref_start
    aligned_query, aligned_ref = "", ""

    while i < len(query) and j < len(ref):
        if query[i] == ref[j]:
            alignment_score += match_score
            aligned_query += query[i]
            aligned_ref += ref[j]
        else:
            mismatches += 1
            alignment_score += mismatch_penalty
            aligned_query += query[i]
            aligned_ref += ref[j]
        i += 1
        j += 1
        if mismatches > mismatch_threshold:
            break

    return aligned_query, aligned_ref, alignment_score, mismatches


def calculate_md_tag(aligned_query, aligned_ref):
    md_tag = ""
    matches = 0
    for q, r in zip(aligned_query, aligned_ref):
        if q == r:
            matches += 1
        else:
            if matches > 0:
                md_tag += str(matches)
            md_tag += r
            matches = 0
    if matches > 0:
        md_tag += str(matches)
    return md_tag


def write_sam_header(ref, output_file):
    with open(output_file, 'w') as f:
        f.write("@HD\tVN:1.0\tSO:unsorted\n")
        for i in range(len(ref.sequence_ids)):
            f.write(f"@SQ\tSN:{ref.sequence_ids[i].replace('>', '')}\tLN:{ref.indexes[ref.sequence_ids[i]]}\n")


def write_sam_record(output_file, query_id, flag, ref_id, ref_start, mapq, cigar, rnext, pnext, tlen, seq, qual, nm, md,
                     as_score, xs):
    with open(output_file, 'a') as f:
        f.write(
            f"{query_id}\t{flag}\t{ref_id}\t{ref_start + 1}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\tNM:i:{nm}\tMD:Z:{md}\tAS:i:{as_score}\tXS:i:{xs}\n")


def main():
    args = parse_arguments()

    start_time = time.time()
    ref = referenceGenome(args.reference_path)
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Indexing reference genome took {execution_time} seconds to execute.")

    start_time = time.time()
    reads = obtain_reads(args.read_path)
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Obtaining reads took {execution_time} seconds to execute.")

    count = 0
    write_sam_header(ref, args.output_file)
    start_time = time.time()
    for query_id in reads:
        count += 1
        query_sequence, query_quality = reads[query_id]
        seeds = seed_and_extend_bwt(query_sequence, ref, args.k, args.mismatch_threshold)
        best_alignment = None
        best_score = float('-inf')
        for seed in seeds:
            query_start, ref_start, length, mismatches = seed

            # Perform full alignment starting from the seeds
            aligned_query, aligned_ref, score, total_mismatches = align(
                query_sequence, ref.sequence, query_start, ref_start, length,
                args.mismatch_penalty, args.gap_penalty, args.match_score, args.mismatch_threshold)

            if score > best_score:
                best_score = score
                md_tag = calculate_md_tag(aligned_query, aligned_ref)
                chrom = ref.chromosome_of_origin(query_start)
                best_alignment = (
                    query_id, 0, chrom, ref_start, 60, f"{length}M", "*", 0, 0,
                    query_sequence, ''.join(chr(q + 33) for q in query_quality), total_mismatches, md_tag,
                    score, 0
                )

        if best_alignment:
            write_sam_record(args.output_file, *best_alignment)
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Seeding and alignment took {execution_time} seconds to execute.")
    print("\n")


if __name__ == "__main__":
    main()
