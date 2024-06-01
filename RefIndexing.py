import numpy as np
import pandas as pd
import os


# think of reference as the hg19 provided in
class referenceGenome:
    def __init__(self, path):
        self.path = path
        self.sequence_id, self.sequence = self.read_reference_genome(self.path)
    def read_reference_genome(self, path:str)-> (str,str):
        with open(path, "r") as file:
            lines = file.readlines()
        sequence_id = lines[0].strip().lstrip('>').rstrip()
        sequence = ''.join(line.strip() for line in lines[1:])
        return sequence_id, sequence
ref =referenceGenome("../data/hg19.fa")
print(ref.sequence_id)
print(ref.sequence[0:30])

