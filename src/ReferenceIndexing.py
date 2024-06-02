import os
from collections import Counter, defaultdict

import numpy as np


class referenceGenome:
    def __init__(self, path):
        self.path = path
        reads = self.read_file()
        self.sequence_ids = list(reads.keys())
        self.sequence = ''.join(reads.values())
        self.indexes = self.data_in()

        # self.sequence_id, self.sequence = self.read_reference_genome(self.path)

        # creating index files + overall index of the sequences provided in the reference
        self.suffix_array = self.suffix_array_construction(self.sequence)
        self.bwt = self.bwt_from_suffix_array(self.sequence, self.suffix_array)
        self.total_counts, self.occ_counts_before = self.build_fm_index(self.bwt)

    # reading fasta file
    def read_file(self):
        samples = {}
        header = None
        with open(self.path, 'r') as file:
            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    header = line
                    samples[header] = ""
                else:
                    samples[header] += line.strip()
        return samples

    def data_seq(self, samples):
        output_file = "./reference.seq"
        out = ""
        if os.path.exists(output_file):
            option = "w"
        else:
            option = "x"
        with open(output_file, option) as file:
            sequences = list(samples.values())
            for i in range(len(sequences)):
                file.write(sequences[i])
                out += sequences[i]
                if i != len(sequences) - 1:
                    file.write("@")
                    out += "@"
        return out

    # indexing file for later use 

    def get_index(self, key, samples, out_seq):
        sequences = out_seq.split("@")
        array_index = sequences.index(samples[key])
        count = 0
        for i in range(array_index):
            count += len(sequences[i]) + 1  # one for @
        return count

    def data_in(self):
        indexes = {}
        samples = self.read_file()
        out_seq = self.data_seq(samples)
        output_file = "./reference.in"
        if os.path.exists(output_file):
            option = "w"
        else:
            option = "x"
        with open(output_file, option) as file:
            for key in samples.keys():
                index = self.get_index(key, samples, out_seq)
                indexes[key] = index
                file.write(f"{key}\t{index}\n")
        return indexes

    def binarySearch(self, i, j, value):
        index_array = list(self.indexes.values())
        while i <= j:
            if i == len(index_array) and j == len(index_array):
                return len(index_array) - 1
            m = (i + j) // 2
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

    def chromosome_of_origin(self, index):
        value_index = self.binarySearch(0, len(list(self.indexes.values())), index)
        if value_index is not None:
            return list(self.indexes.keys())[value_index].split(" ")[0].replace(">", "")
        else:
            return "unknown"

    def suffix_array_construction(self, s: str) -> np.ndarray:
        """Constructs the suffix array of a given string s using SA-IS algorithm."""
        s = s + "$"
        n = len(s)
        suffix_array = np.zeros(n, dtype=int)

        def sort_cyclic_shifts():
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

        sort_cyclic_shifts()
        return np.array(suffix_array)

    # Function to create BWT from suffix array
    def bwt_from_suffix_array(self, s: str, suffix_array: np.ndarray) -> np.ndarray:
        s = s + "$"
        bwt = ''.join(s[i - 1] if i > 0 else s[-1] for i in suffix_array)
        return bwt

    # Function to build the FM-index (C array and O table)
    def build_fm_index(self, bwt: np.ndarray):
        # C array
        counts = Counter(bwt)
        total_counts = dict()
        sum_counts = 0
        for char in sorted(counts.keys()):
            total_counts[char] = sum_counts
            sum_counts += counts[char]

        # O table
        occ_counts_before = defaultdict(lambda: [0] * (len(bwt) + 1))
        for i in range(len(bwt)):
            char = bwt[i]
            for c in occ_counts_before:
                occ_counts_before[c][i + 1] = occ_counts_before[c][i]
            occ_counts_before[char][i + 1] += 1

        return total_counts, occ_counts_before

    # Function for backward search using FM-index
    def backward_search(self, pattern):
        l, r = 0, len(self.bwt) - 1
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