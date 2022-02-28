import os
from collections import OrderedDict

__all__ = ["GFFParser"]

class GFFParser:
    def __init__(self, gff):
        self.gff = os.path.basename(gff)
        self.path = os.path.realpath(gff)
        self.gff_name = "_".join(self.gff.strip(".gff").split("_")[:3])

    def parse_gff(self):
        with open(self.gff, "r") as fd:
            for row in fd.readlines():
                if row == "##FASTA\n":
                    break
                elif "#" in row:
                    continue
                else:
                    row = row.split("\t")
                    d = row[-1].split(";")
                    d = OrderedDict([(s.split("=")[0], s.split("=")[1]) for s in d])
                    row[-1] = d
                    yield row
    
