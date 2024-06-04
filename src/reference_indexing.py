import os
from collections import Counter, defaultdict

import numpy as np


class ReferenceGenome:
    """Class for the reference genome object.

    Attributes:
        path (str): The path to the reference genome file.
        sequence_ids (list[str]): A list of the sequence IDs.
        sequence (str): The concatenated sequences.
        indexes (dict[str, int]): A dictionary of the sequence IDs and indexes.
        suffix_array (np.ndarray): The suffix array of the reference genome.
        bwt (str): The Burrows-Wheeler Transform (BWT) of the reference genome.
        total_counts (dict[str, int]): The total counts of characters in the BWT.
        occ_counts_before (defaultdict[str, list[int]]): The occurrence counts before a character in the BWT.

    Methods:
        read_file: Reads the reference genome fasta file.
        data_seq: Writes the sequences to a file.
        get_index: Gets the index of the sequence in the concatenated sequences.
        data_in: Writes the indexes of the sequences to a file.
        binarySearch: Performs a binary search on the indexes of the sequences.
        chromosome_of_origin: Gets the chromosome of origin of a sequence given its index.
        sort_cyclic_shifts: Sorts the cyclic shifts of a string s using SA-IS algorithm.
        suffix_array_construction: Constructs the suffix array of a given string s using SA-IS algorithm.
        bwt_from_suffix_array: Creates the Burrows-Wheeler Transform (BWT) of a given string s using the suffix array.
        build_fm_index: Builds the FM-index (C array and O table).
        backward_search: Performs backward search using FM-index.
    """

    def __init__(self, path):
        """Constructor for the ReferenceGenome class.

        Args:
            path (str): The path to the reference genome file.
        """
        self.path: str = path
        reads: dict[str, str] = self.read_file()
        self.sequence_ids: list[str] = list(reads.keys())
        self.sequence: str = ''.join(reads.values())
        self.indexes: dict[str, int] = self.data_in()

        # self.sequence_id, self.sequence = self.read_reference_genome(self.path)

        # creating index files + overall index of the sequences provided in the reference
        self.suffix_array: np.ndarray = self.suffix_array_construction(self.sequence)
        self.bwt: str = self.bwt_from_suffix_array(self.sequence, self.suffix_array)

        # building FM-index
        self.total_counts: dict[str, int]
        self.occ_counts_before: defaultdict[str, list[int]]
        self.total_counts, self.occ_counts_before = self.build_fm_index(self.bwt)

    # reading fasta file
    def read_file(self) -> dict[str, str]:
        """Reads the reference genome fasta file. From Wikipedia: A sequence 
        begins with a greater-than character (">") followed by a description 
        of the sequence (all in a single line). The lines immediately following
        the description line are the sequence representation, with one letter 
        per amino acid or nucleic acid, and are typically no more than 80 
        characters in length.

        Returns:
            dict[str, str]: A dictionary of the sequence IDs and sequences.
        """
        samples: dict[str, str] = {}
        header: str = ""
        with open(self.path, 'r') as file:
            for line in file:
                line: str = line.strip()
                if line.startswith(">"):
                    header = line
                    samples[header] = ""
                else:
                    samples[header] += line.strip()
        return samples

    def data_seq(self, samples: dict[str, str]) -> str:
        """Writes the sequences to a file.

        Args:
            samples (dict[str, str]): A dictionary of the sequence IDs and sequences.

        Returns:
            str: The concatenated sequences.
        """
        output_file: str = "./reference.seq"
        out: str = ""
        option: str
        if os.path.exists(output_file):
            option = "w"
        else:
            option = "x"
        with open(output_file, option) as file:
            sequences: list[str] = list(samples.values())
            for i in range(len(sequences)):
                file.write(sequences[i])
                out += sequences[i]
                if i != len(sequences) - 1:
                    file.write("@")
                    out += "@"
        return out

    def get_index(self, key: str, samples: dict[str, str], out_seq: str) -> int:
        """Gets the index of the sequence in the concatenated sequences.

        Args:
            key (str): Sequence ID of sequence being indexed.
            samples (dict[str, str]): A dictionary of the sequence IDs and sequences.
            out_seq (str): The concatenated sequences.

        Returns:
            int: The index of the sequence in the concatenated sequences.
        """
        sequences: list[str] = out_seq.split("@")
        array_index: int = sequences.index(samples[key])
        count: int = 0
        for i in range(array_index):
            count += len(sequences[i]) + 1  # one for @
        return count

    def data_in(self) -> dict[str, int]:
        """Writes the indexes of the sequences to a file.

        Returns:
            dict[str, int]: A dictionary of the sequence IDs and indexes.
        """
        indexes: dict[str, int] = {}
        samples: dict[str, str] = self.read_file()
        out_seq: str = self.data_seq(samples)
        output_file: str = "./reference.in"
        option: str
        if os.path.exists(output_file):
            option = "w"
        else:
            option = "x"
        with open(output_file, option) as file:
            for key in samples.keys():
                index: int = self.get_index(key, samples, out_seq)
                indexes[key] = index
                file.write(f"{key}\t{index}\n")
        return indexes

    def binarySearch(self, i: int, j: int, value: int) -> int | None:
        """Performs a binary search on the indexes of the sequences.

        Args:
            i (int): Lower bound of the search.
            j (int): Upper bound of the search.
            value (int): The value to search for.

        Returns:
            int | None: The index of the sequence in the concatenated sequences,
            if it exists.
        """
        index_array: list[int] = list(self.indexes.values())
        while i <= j:
            if i == len(index_array) and j == len(index_array):
                return len(index_array) - 1
            m: int = (i + j) // 2
            if m == len(index_array):
                m = len(index_array) - 1
            if index_array[m] == value:
                return m
            if index_array[m - 1] <= value and value < index_array[m]:
                return m - 1
            elif index_array[m] > value:
                j = m - 1
            else:
                i = m + 1
        return None

    def chromosome_of_origin(self, index: int) -> str:
        """Gets the chromosome of origin of a sequence given its index.

        Args:
            index (int): The index of the sequence in the concatenated sequences.

        Returns:
            str: The chromosome of origin of the sequence, or "unknown".
        """
        value_index: int | None = self.binarySearch(0, len(list(self.indexes.values())), index)
        if value_index is not None:
            return list(self.indexes.keys())[value_index].split(" ")[0].replace(">", "")
        else:
            return "unknown"

    def sort_cyclic_shifts(self, s: str, n: int, suffix_array: np.ndarray) -> None:
        """Sorts the cyclic shifts of a string s using SA-IS algorithm.

        Args:
            s (str): The input string.
            n (int): The length of the input string.
            suffix_array (np.ndarray): The suffix array of the input string.
        """
        alphabet = 256
        p = np.zeros(n, dtype=int)
        c = np.zeros(n, dtype=int)
        cnt = np.zeros(max(alphabet, n), dtype=int)

        for i in range(n):
            cnt[ord(s[i])] += 1
        for i in range(1, alphabet):
            cnt[i] += cnt[i - 1]
        for i in range(n - 1, -1, -1):
            cnt[ord(s[i])] -= 1
            p[cnt[ord(s[i])]] = i
        classes = 1
        c[p[0]] = 0
        for i in range(1, n):
            if s[p[i]] != s[p[i - 1]]:
                classes += 1
            c[p[i]] = classes - 1

        pn = np.zeros(n, dtype=int)
        cn = np.zeros(n, dtype=int)
        h = 0
        while (1 << h) < n:
            for i in range(n):
                pn[i] = p[i] - (1 << h)
                if pn[i] < 0:
                    pn[i] += n
            cnt = np.zeros(classes, dtype=int)
            for i in range(n):
                cnt[c[pn[i]]] += 1
            for i in range(1, classes):
                cnt[i] += cnt[i - 1]
            for i in range(n - 1, -1, -1):
                cnt[c[pn[i]]] -= 1
                p[cnt[c[pn[i]]]] = pn[i]
            classes = 1
            cn[p[0]] = 0
            for i in range(1, n):
                curr = (c[p[i]], c[(p[i] + (1 << h)) % n])
                prev = (c[p[i - 1]], c[(p[i - 1] + (1 << h)) % n])
                if curr != prev:
                    classes += 1
                cn[p[i]] = classes - 1
            c, cn = cn, c
            h += 1

        for i in range(n):
            suffix_array[i] = p[i]

    def suffix_array_construction(self, s: str) -> np.ndarray:
        """Constructs the suffix array of a given string s using SA-IS algorithm.
        
        Args:
            s (str): The input string.

        Returns:
            np.ndarray: The suffix array of the input string.
        """
        s: str = s + "$"
        n: int = len(s)
        suffix_array: np.ndarray = np.zeros(n, dtype=int)
        self.sort_cyclic_shifts(s, n, suffix_array)
        return np.array(suffix_array)

    def bwt_from_suffix_array(self, s: str, suffix_array: np.ndarray) -> str:
        """Creates the Burrows-Wheeler Transform (BWT) of a given string s using
        the suffix array.

        Args:
            s (str): The input string.
            suffix_array (np.ndarray): The suffix array of the input string.

        Returns:
            str: The BWT of the input string.
        """
        s: str = s + "$"
        bwt: str = ''.join(s[i - 1] if i > 0 else s[-1] for i in suffix_array)
        return bwt

    def build_fm_index(self, bwt: np.ndarray) -> tuple[dict[str, int], 
                                                       defaultdict[str, list[int]]]:
        """Builds the FM-index (C array and O table).

        Args:
            bwt (np.ndarray): The Burrows-Wheeler Transform (BWT) of the reference genome.

        Returns:
            tuple[dict[str, int], defaultdict[str, list[int]]]: A tuple containing the total counts
        """
        # C array
        counts: Counter[str] = Counter(bwt)
        total_counts: dict[str, int] = dict()
        sum_counts: int = 0
        for char in sorted(counts.keys()):
            total_counts[char] = sum_counts
            sum_counts += counts[char]

        # O table
        occ_counts_before: defaultdict[str, list[int]]
        occ_counts_before = defaultdict(lambda: [0] * (len(bwt) + 1))

        for i in range(len(bwt)):
            char: str = bwt[i]
            for c in occ_counts_before:
                occ_counts_before[c][i + 1] = occ_counts_before[c][i]
            occ_counts_before[char][i + 1] += 1

        return total_counts, occ_counts_before

    # Function for backward search using FM-index
    def backward_search(self, pattern: str) -> np.ndarray:
        """Performs backward search using FM-index.

        Args:
            pattern (str): The pattern to search for.

        Returns:
            np.ndarray: The positions of the pattern in the reference genome.
        """
        l: int = 0
        r: int = len(self.bwt) - 1
        for char in reversed(pattern):
            if char == 'N':
                continue  # Skip 'N' characters
            try:
                l = self.total_counts[char] + self.occ_counts_before[char][l]
                r = self.total_counts[char] + self.occ_counts_before[char][r + 1] - 1
            except KeyError:
                return []  # Return empty list if char is not found in total_counts
            if l > r:
                return []
        return self.suffix_array[l:r + 1]
