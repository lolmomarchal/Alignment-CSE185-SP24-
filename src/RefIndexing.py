from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import os
import array


# think of reference as the hg19 provided in
class referenceGenome:
    def __init__(self, path):
        self.path = path
        self.sequence_id, self.sequence = self.read_reference_genome(self.path)
        self.array = self.suffix_array_construction(self.sequence)
        self.bwt = self.bwt_from_suffix_array(self.sequence, self.array)
        self.total_counts, self.occ_counts_before = self.build_fm_index(self.bwt)

    def read_reference_genome(self, path:str)-> (str,str):
        with open(path, "r") as file:
            lines = file.readlines()
        sequence_id = lines[0].strip().lstrip('>').rstrip()
        sequence = ''.join(line.strip() for line in lines[1:])
        return sequence_id, self.filter_unknown_regions(sequence)
    def filter_unknown_regions(self, sequence:str)->str:
        return sequence.replace("N", "")
    def suffix_array_construction(self,s):
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
        return suffix_array

    # Function to create BWT from suffix array
    def bwt_from_suffix_array(self,s, suffix_array):
        s = s + "$"
        bwt = ''.join(s[i - 1] if i > 0 else s[-1] for i in suffix_array)
        return bwt

    # Function to build the FM-index (C array and O table)
    def build_fm_index(self,bwt):
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
    def backward_search(self,pattern, bwt, suffix_array, total_counts, occ_counts_before):
        l, r = 0, len(bwt) - 1
        for char in reversed(pattern):
            l = total_counts[char] + occ_counts_before[char][l]
            r = total_counts[char] + occ_counts_before[char][r + 1] - 1
            if l > r:
                return []
        return suffix_array[l:r + 1]




ref =referenceGenome("../data/hg19.fa")

output = ref.backward_search("GTACGAGCTCAACGACGGTAAAGAGGAT", ref.bwt, ref.array,ref.total_counts, ref.occ_counts_before)
print(output)
