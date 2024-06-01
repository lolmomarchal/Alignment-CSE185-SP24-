# outside imports
from Bio import SeqIO
import time
# inside imports
from ReferenceIndexing import referenceGenome as referenceGenome


def main():
    start_time = time.time()
    ref = referenceGenome("../data/hg19.fa")
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Indexing reference genome took {execution_time} seconds to execute.")


if __name__ == "__main__":
    main()