import os
import numpy
import json
import subprocess
import csv

from collections import OrderedDict

from genomictools.wrapper import timer
from genomictools.vcf.VCFParser import VCFParser
from genomictools.genbank.Genome import Genome


__all__ = ["SNPMatrix"]

class SNPMatrix:
    @timer
    def __init__(self, vcf=None, ref_genome=None, matrix=None, ref_name=None, strains_index=None):
        if matrix and ref_name and strains_index:
            self.ref_name = ref_name
            if isinstance(matrix, str):
                self.matrix = self.load_matrix(matrix)
            elif isinstance(matrix, numpy.ndarray):
                self.matrix = matrix
            else:
                exit("matrix format is not supported...")
            self.index = strains_index
        elif vcf and ref_genome:
            self.ref_genome = Genome(ref_genome)
            self.ref_name = self.ref_genome.name
            self.vcfparser = VCFParser(vcf)
            self.matrix = numpy.array(self.build_matrix())
            self.index = list(self.strains_index())
        else:
            print("Give a VCF and a genbank file for ref OR a matrix already build with ref_name and strain index [(i, strain), (i, strain)]")

    def strain_index(self):
        #call = subprocess.Popen(["vcfsamplenames", self.combined_vcf], stdout=subprocess.PIPE)
        #data = call.communicate()[0].split(b"\n")
        #data = [b.decode("utf8") for b in data if b.decode("utf8")]
        for i, strain in enumerate(self.vcfparser.samples):
            yield (i, strain)

    def link_snp_to_locustag(self, pos):
        for locustag, feature in self.ref_genome:
            if feature.start <= pos <= feature.end:
                return locustag
        return "Intergenic"

    def build_matrix(self):
        matrix = []
        #plain_motif = ".:.:.:.:.:.:.:."
        #empty_motif = "."
        #vcftools_motif = "./."
        for row in self.vcfparser.parse_snp():
            line = [row.pos, row.ref]
            for sample, format in row.samples.items():
                if format == "NA":
                    line.append("N")
                elif format["GT"][0] == "0":
                    line.append(row.ref)
                else:
                    line.append(row.alt)
            #for i in range(9, len(row.row)):
            #    if row.row[i] == plain_motif or row.row[i] == empty_motif or row.row[i].split(":")[0] == vcftools_motif:
            #        line.append("N")
            #    elif row.row[i].split(":")[0] == "0":
            #        line.append(row.ref)
            #    else:
            #        line.append(row.alt)
            annotation = self.link_snp_to_locustag(int(row.pos))
            line.append(annotation)
            matrix.append(line)
        return matrix

    def build_snp_concatenates_from_matrix(self):
        cat = ">" + self.ref_name + "\n"
        cat += "".join(self.matrix[:,1])
        cat += "\n"
        for (pos, strain) in self.index:
            cat += ">" + strain + "\n"
            cat += "".join(self.matrix[:,pos + 2])
            cat += "\n"
        return cat

    def save_matrix(self, output):
        with open(output, 'w') as fd:
            json.dump(self.matrix.tolist(), fd, indent=4)

    def load_matrix(self, json_matrix):
        with open(json_matrix, 'r') as fd:
            matrix = json.load(fd)
        return numpy.array(matrix)

    def save_concatenate(self, output):
        with open(output, "w") as fd:
            fd.write(self.build_snp_concatenates_from_matrix())

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        rows = []

        head = ["Position", self.ref_name]
        for i, strain in self.index:
            head.append(strain)
        head.append("Annotation")
        head = "\t".join(head[0:3] + ["..."] + head[-3:-1] + [head[-1]])
        rows.append(head)

        for n in range(0, 10, 1):
            row = self.matrix[n,:]
            rows.append("\t".join(list(row[0:3]) + ["..."] + list(row[-3:-1]) + [row[-1]]))
        rows.append("...")

        for n in range(len(self.matrix) - 10, len(self.matrix), 1):
            row = self.matrix[n,:]
            rows.append("\t".join(list(row[0:3]) + ["..."] + list(row[-3:-1]) + [row[-1]]))
        rows = "\n".join(rows)
        
        print(rows)

        return f"""
        Reference genome: {self.ref_name}
        Number of strain: {len(self.index) + 1}
        Matrix length: {len(self.matrix)}
        """
