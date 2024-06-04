import argparse
import time

from Bio import SeqIO

from reference_indexing import ReferenceGenome

# Parameters for seed and extend
DEFAULT_K: int = 10
DEFAULT_MISMATCH_THRESHOLD: int = 3
DEFAULT_GAP_PENALTY: int = -2
DEFAULT_MISMATCH_PENALTY: int = -1
DEFAULT_MATCH_SCORE: int = 2


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Seed and Extend BWT Aligner")
    parser.add_argument('-k', '--k', type=int, default=DEFAULT_K, 
                        help=f'Length of k-mer for seeding (Default: {DEFAULT_K})')
    parser.add_argument('-mt', '--mismatch_threshold', type=int, default=DEFAULT_MISMATCH_THRESHOLD,
                        help=f'Threshold for mismatches (Default: {DEFAULT_MISMATCH_THRESHOLD})')
    parser.add_argument('-d', '--gap_penalty', type=int, default=DEFAULT_GAP_PENALTY, 
                        help=f'Penalty for gaps (Default: {DEFAULT_GAP_PENALTY})')
    parser.add_argument('-s', '--mismatch_penalty', type=int, default=DEFAULT_MISMATCH_PENALTY,
                        help=f'Penalty for mismatches (Default: {DEFAULT_MISMATCH_PENALTY})')
    parser.add_argument('-m', '--match_score', type=int, default=DEFAULT_MATCH_SCORE, 
                        help=f'Score for matches (Default: {DEFAULT_MATCH_SCORE})')
    parser.add_argument('-r', '--reference_path', type=str, required=True, help='Path to reference genome file')
    parser.add_argument('-i', '--read_path', type=str, required=True, help='Path to reads file')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to output SAM file')
    return parser.parse_args()


def seed_and_extend_bwt(query: str, reference: ReferenceGenome, 
                        k: int, mismatch_threshold: int
                        ) -> list[tuple[int, int, int, int]]:
    """Seed and extend using the Burrows-Wheeler Transform (BWT) algorithm.

    Args:
        query (str): The query sequence.
        reference (ReferenceGenome): The reference genome object.
        k (int): The length of the k-mer for seeding.
        mismatch_threshold (int): The threshold for mismatches.

    Returns:
        list[tuple[int, int, int, int]]: A list of tuples containing the 
        seed positions.
    """
    seeds: list[tuple[int, int, int, int]] = []
    for i in range(len(query) - k + 1):
        kmer: str = query[i:i + k]
        positions = reference.backward_search(kmer)
        for pos in positions:
            match_length: int
            mismatches: int
            match_length, mismatches = extend_match(query, reference.sequence, i, pos, k, mismatch_threshold)
            if match_length > k:
                seeds.append((i, pos, match_length, mismatches))
    return seeds


def extend_match(query: str, reference: ReferenceGenome, query_start: int,
                 ref_start: int, seed_length: int, mismatch_threshold: int
                ) -> tuple[int, int]:
    """Extend the match between the query and reference sequences.

    Args:
        query (str): The query sequence.
        reference (ReferenceGenome): The reference genome object.
        query_start (int): Index of the query sequence to start extending.
        ref_start (int): Index of the reference sequence to start extending.
        seed_length (int): The length of the seed.
        mismatch_threshold (int): The threshold for mismatches.

    Returns:
        tuple[int, int]: A tuple containing the match length and the number of mismatches.
    """
    match_length: int = seed_length
    mismatches: int = 0
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


def obtain_reads(path_reads: str) -> dict[str, tuple[str, list[int]]]:
    """Obtain reads from a FASTQ file.

    Args:
        path_reads (str): The path to the FASTQ file containing the reads.

    Returns:
        dict[str, tuple[str, list[int]]]: A dictionary containing the reads.
    """
    sequences: dict[str, tuple[str, list[int]]] = {}
    for record in SeqIO.parse(path_reads, "fastq"):
        record: SeqIO.SeqRecord
        sequences[record.id] = (str(record.seq), record.letter_annotations["phred_quality"])
    return sequences


def align(query: str, ref: ReferenceGenome, query_start: int, ref_start: int,
          match_length: int, mismatch_penalty: int, gap_penalty: int, 
          match_score: int, mismatch_threshold: int) -> tuple[str, str, int, int]:
    """Align the query and reference sequences.

    Args:
        query (str): The query sequence.
        ref (ReferenceGenome): The reference genome object.
        query_start (int): Index of the query sequence to start aligning.
        ref_start (int): Index of the reference sequence to start aligning.
        match_length (int): The length of the match.
        mismatch_penalty (int): Penalty for mismatches.
        gap_penalty (int): Penalty for gaps.
        match_score (int): Score for matches.
        mismatch_threshold (int): The threshold for mismatches.

    Returns:
        tuple[str, str, int, int]: A tuple containing the aligned query,
        aligned reference, alignment score, and the number of mismatches.
    """
    alignment_score: int = 0
    mismatches: int = 0
    i, j = query_start, ref_start
    aligned_query: str
    aligned_ref: str
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


def calculate_md_tag(aligned_query: str, aligned_ref: str) -> str:
    """Calculate the MD tag for the SAM file, the string encoding
    the mismatched and deleted reference bases.

    Args:
        aligned_query (str): The aligned query sequence.
        aligned_ref (str): The aligned reference sequence.

    Returns:
        str: The MD tag for the SAM file.
    """
    md_tag: str = ""
    matches: int = 0
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


def write_sam_header(ref: ReferenceGenome, output_file: str) -> None:
    """Write the headers to the output SAM file.

    Args:
        ref (ReferenceGenome): The reference genome object.
        output_file (str): The path to the output SAM file.
    """
    with open(output_file, 'w') as f:
        f.write("@HD\tVN:1.0\tSO:unsorted\n")
        for i in range(len(ref.sequence_ids)):
            f.write(f"@SQ\tSN:{ref.sequence_ids[i].replace('>', '')}\tLN:{ref.indexes[ref.sequence_ids[i]]}\n")


def write_sam_record(output_file: str, query_id: str, flag: int, ref_id: str, ref_start: int, 
                     mapq: int, cigar: str, rnext: str, pnext: int, tlen: int, seq: str, 
                     qual: str, nm: int, md: str, as_score: int, xs: int) -> None:
    """Write a SAM record to the output file.

    Args:
        output_file (str): The path to the output SAM file.
        query_id (str): SAM record query ID.
        flag (int): SAM record flag.
        ref_id (str): SAM record reference ID.
        ref_start (int): SAM record reference start.
        mapq (int): SAM record mapping quality.
        cigar (str): SAM record CIGAR string.
        rnext (str): SAM record rnext.
        pnext (int): SAM record pnext.
        tlen (int): SAM record tlen.
        seq (str): SAM record sequence.
        qual (str): SAM record quality.
        nm (int): Number of mismatches.
        md (str): MD tag.
        as_score (int): Alignment score.
        xs (int): Alignment score for secondary alignments.
    """
    with open(output_file, 'a') as f:
        f.write(
            "\t".join(str(x) for x in [query_id, flag, ref_id, ref_start + 1, 
            mapq, cigar, rnext, pnext, tlen, seq, qual, f"NM:i:{nm}", 
            f"MD:Z:{md}", f"AS:i:{as_score}", f"XS:i:{xs}"]) + "\n"
        )


def main() -> None:
    args: argparse.Namespace = parse_arguments()

    start_time: float = time.time()
    ref: ReferenceGenome = ReferenceGenome(args.reference_path)
    end_time: float = time.time()
    execution_time: float = end_time - start_time
    print(f"Indexing reference genome took {execution_time} seconds to execute.")

    start_time: float = time.time()
    reads: dict[str, tuple[str, list[int]]] = obtain_reads(args.read_path)
    end_time: float = time.time()
    execution_time: float = end_time - start_time
    print(f"Obtaining reads took {execution_time} seconds to execute.")

    count: int = 0
    write_sam_header(ref, args.output_file)
    start_time: float = time.time()
    for query_id in reads:
        
        query_id: str
        query_sequence: str
        query_quality: tuple[str, list[int]]
        seeds: list[tuple[int, int, int, int]]

        count += 1
        query_sequence, query_quality = reads[query_id]
        seeds = seed_and_extend_bwt(query_sequence, ref, args.k, args.mismatch_threshold)
        best_alignment: tuple | None = None
        best_score: float = float('-inf')
        for seed in seeds:

            query_start: int
            ref_start: int
            length: int
            mismatches: int

            query_start, ref_start, length, mismatches = seed

            # Perform full alignment starting from the seeds
            aligned_query: str
            aligned_ref: str
            score: int
            total_mismatches: int

            aligned_query, aligned_ref, score, total_mismatches = align(
                query_sequence, ref.sequence, query_start, ref_start, length,
                args.mismatch_penalty, args.gap_penalty, args.match_score, args.mismatch_threshold)

            if score > best_score:
                best_score = score
                md_tag: str = calculate_md_tag(aligned_query, aligned_ref)
                chrom: str = ref.chromosome_of_origin(query_start)
                best_alignment = (
                    query_id, 0, chrom, ref_start, 60, f"{length}M", "*", 0, 0,
                    query_sequence, ''.join(chr(q + 33) for q in query_quality), total_mismatches, md_tag,
                    score, 0
                )

        if best_alignment:
            write_sam_record(args.output_file, *best_alignment)

    end_time: float = time.time()
    execution_time: float = end_time - start_time
    print(f"Seeding and alignment took {execution_time} seconds to execute.")
    print("\n")


if __name__ == "__main__":
    main()
