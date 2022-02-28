import os
import csv
import numpy
import pandas
import re
from io import StringIO

from genomitools.wrapper import timer, pprint

__all__ = ["VCFParser", "VCFParserOrigin", "VCFParserNew"]

def VCFParser(vcf, method="origin"):
    if method == "origin":
        return VCFParserOrigin(vcf)
    else:
        return VCFParserNew(vcf)


from collections import OrderedDict
import numpy
class VCFParserNew:
    #@profile
    #@timer
    def __init__(self, vcf):
        self.vcf = os.path.realpath(vcf)
        self.header = None
        self.samples = None
        self.content = None
        self.samples_data = None
        self.data = None
        self.read_vcf()
        self.build()

    #@profile
    def read_vcf(self):
        with open(self.vcf, "r") as fd:
            for header in fd:
                if "#CHROM" in header:
                    self.header = header.strip().split("\t")
                    self.samples = {sample: i for i, sample in enumerate(self.header[9:])}
                    break
            fd.seek(0) # readlines parcours le fichier complètement donc on doit se remettre au début
            self.content = numpy.array([row.strip().split("\t") for row in fd if row[0] != "#"])
            self.samples_data = self.content[:, 9:]
            print(self.content.shape)

    def build(self):
        self.data = OrderedDict()
        self.data["CHROM"] = self.content[:, 0]
        self.data["POS"] = numpy.array(self.content[:, 1], dtype=numpy.uint32)
        self.data["REF"] = self.content[:, 3]
        self.data["ALT"] = self.content[:, 4]
        self.data["QUAL"] = numpy.array(self.content[:, 5], dtype=numpy.float32)
        self.data["INFO"] = self.content[:, 7]

    def buildGenotypeArray(self):
        return GenotypeArray(self)

# From scikit-alleles
class GenotypeArray:
    #@timer
    def __init__(self, parser):
        self.samplesIndex = {sample: index for index, sample in enumerate(parser.samples)}
        self.formatIndex = {f: index for index, f in enumerate(['GT', 'DP', 'AD', 'RO', 'QR', 'AO', 'QA', 'GL'])}
        self.samplesFormat = numpy.apply_along_axis(self.splitFormat, 1, parser.samples_data)
        self.genotype = pandas.DataFrame(data=self.samplesFormat[:, :, 0], columns=self.samplesIndex.keys()).apply(self.convertGenotype)
        self.depth = pandas.DataFrame(data=self.samplesFormat[:, :, 1], columns=self.samplesIndex.keys()).apply(self.convertDepth)

    def splitFormat(self, col, *args, **kwargs):
        return [c.split(":") for c in col]

    def convertGenotype(self, col, *args, **kwargs):
        return [-1 if i == "./." else int(i) for i in col]

    def convertDepth(self, col, *args, **kwargs):
        return [numpy.nan if i == "." else int(i) for i in col]

class VCFParserOrigin:
    def __init__(self, vcf):
        self.vcf = os.path.realpath(vcf)
        self.metadata = dict()
        self.contig = dict()
        self.infofields = dict()
        self.formatfields = dict()
        self.filterfields = dict()
        self.convert_type = {"Integer": int, "Float": float, "String": str, "Character": str, "Flag": True}
        self.header = []
        self._parse_header()
        if self.header:
            self.samples = self.header[9:] # subset of header
        else:
            self.samples = []

    def _get_header_line(self):
        with open(self.vcf, "r") as fd:
            return [row.strip() for row in fd if row[0] == "#"]

    def _search_field(self, row):
        results = []
        # From vcflib : https://github.com/vcflib/vcflib/blob/40dbb399b5d25ae694e15755724475b274d1b8fe/scripts/vcf2sqlite.py
        i = re.search('ID=(.*?),', row)
        n = re.search('Number=(.*?),', row)
        t = re.search('Type=(.*?),', row)
        d = re.search('Description="(.*?)"', row)
        if i and n and t and d:
            results.append(i.groups()[0])
            results.append(n.groups()[0])
            results.append(t.groups()[0])
            results.append(d.groups()[0])
            return results
        else:
            raise ValueError(f"Error while parsing rows:\n{row}\n")
    
    def _parse_header(self):
        for row in self._get_header_line():
            if row[:2] == "##":
                if row[:8] == "##contig":
                    i = re.search('ID=(.*?),', row)
                    l = re.search('length=(.*?)>', row)
                    if i and l:
                        self.contig["name"] = i.groups()[0]
                        self.contig["length"] = int(l.groups()[0])
                    else:
                        print("Error in contig search")
                    #i = re.search('ID=(.*?),', row)
                    #l = re.search('length=(.*?),', row)
                    #if i and l:
                    #    key = i.groups()[0]
                    #    length = l.groups()[0]
                    #    value = {"length": length}
                    #    self.contig[key] = value
                elif row[:6] == "##INFO":
                    key, number, type, description = self._search_field(row)
                    value = {"Number": number, "Type": type, "Description": description}
                    self.infofields[key] = value
                elif row[:8] == "##FORMAT":
                    key, number, type, description = self._search_field(row)
                    value = {"Number": number, "Type": type, "Description": description}
                    self.formatfields[key] = value
                elif row[:8] == "##FILTER":
                    i = re.search('ID=(.*?),', row)
                    d = re.search('Description="(.*?)"', row)
                    if i and d:
                        key = i.groups()[0]
                        description = d.groups()[0]
                        value = {"Description": description}
                        self.filterfields[key] = value
                else:
                    row = row.strip("#").replace('"', '')
                    key = row.split("=")[0]
                    value = row.split("=")[1]
                    self.metadata[key] = value
            elif row[:6] == "#CHROM":
                row = row.split("\t")
                row[0] = row[0].strip("#")
                self.header = row
                break
            else:
                raise Exception(f"Error while parsing header of vcf file: {self.vcf}\nUnsupported file !")

    def read_vcf(self):
        with open(self.vcf, 'r') as fd:
            lines = [l for l in fd if not l.startswith('##')]
        return pandas.read_csv(StringIO(''.join(lines)), 
                dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str, 
                        'QUAL': float, 'FILTER': str, 'INFO': str, 'FORMAT': str}, 
                sep='\t').rename(columns={'#CHROM': 'CHROM'})

    # Takes on input each element of the column individualy and return the modified element
    #@staticmethod
    #def _apply_to_database(elem):
    #    if elem.split(":")[0] in [".", "./."]:
    #        return "NA"
    #    else:
    #        return elem.split(":")

    def parse_snp(self):
        with open(self.vcf, "r") as fd:
            for row in fd:
                if not row[0] == "#":
                    row = row.strip().split("\t")
                    if len(row[3]) == len(row[4]) == 1 and "TYPE=snp" in row[7]:
                        yield VCFRow(self, row)


class VCFRow:
    @timer
    def __init__(self, parser, row):
        self.parser = parser
        self.row = row
        self.samples = dict() 
        self._set_attr()
        self.set_samples()

    def __getitem__(self, key):
        return self.__getattribute__(key)

    def __repr__(self):
        return f"{self.parsed_row()}"

    def __str__(self):
        return self.__repr__()

    def _set_attr(self):
        for i, col in enumerate(self.row[:9]):
            if i == 7: # INFO field
                pairs = col.split(";")
                result = dict()
                for pair in pairs:
                    key = pair.split("=")[0]
                    number = self.parser.infofields[key]["Number"]
                    t = self.parser.convert_type[self.parser.infofields[key]["Type"]]
                    if number == "0":
                        result[key] = True
                    elif number == "1":
                        result[key] = t(pair.split("=")[1])
                    elif number > "1" or number in ["A", "R", "G", "."]:
                        val = pair.split("=")[1]
                        result[key] = [t(v) for v in val.split(",")]
                        self.__setattr__(self.parser.header[i].lower(), result)
            elif i == 8: # FORMAT
                self.__setattr__(self.parser.header[i].lower(), col.split(":"))
            else:
                self.__setattr__(self.parser.header[i].lower(), col)
        else:
            self.pos = int(self.pos)
            self.qual = float(self.qual)

    def set_samples(self):
        #self.samples = {sample: {k: [float(n) for n in v] for k, v in zip(self.format, [f.split(",") for f in format])} if format else "NA" for sample, format in zip(self.parser.samples, [elem.split(":") if not "." in elem.split(":")[0] else None for elem in self.row[9:]])}
        for sample, format in zip(self.parser.samples, [elem.split(":") if not "." in elem.split(":")[0] else None for elem in self.row[9:]]):
            if format:
                format = [f.split(",") for f in format]
                self.samples[sample] = {k: [float(n) for n in v] for k, v in zip(self.format, format)}
            else:
                self.samples[sample] = "NA"

    def parsed_row(self):
        parsed_row = [col if i < 9 else {self.parser.samples[i - 9]: self.samples[self.parser.samples[i - 9]]} for i, col in enumerate(self.row)]
        return str(parsed_row)

    def depth_recalibration(self):
        DP = [format["DP"][0] for format in self.samples.values() if format != "NA"]
        try:
            info_depth = sum(DP) // len(DP)
            return info_depth
        except ZeroDivisionError:
            print(f"{self.pos} -> {DP}")
            return 0


