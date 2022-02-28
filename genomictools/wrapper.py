"""
This module define some personnal wrapper.
"""
import time
import subprocess
import os

from pprint import PrettyPrinter
pprint = PrettyPrinter(indent=4).pprint
del PrettyPrinter

__all__ = ["pprint", "timer", "blastn", "blastp", "BlastResult", "Alignment", "check"]

class Alignment:
    def __init__(self, fields, row):
        self.fields = fields
        self.row = row
        for n, field in enumerate(row):
            name = self.fields[n].replace(' ', '_')
            try:
                self.__setattr__(name, float(field))
            except Exception:
                self.__setattr__(name, field)

class BlastResult:
    def __init__(self, path):
        self.path = path
        self.result = []
        with open(self.path, "r") as fd:
            self.data = fd.read()
            for n, row in enumerate(self.data.splitlines()):
                if not row[0] == "#":
                    row = row.split('\t')
                    self.result.append(Alignment(self.fields, row))
                else:
                    if n == 0:
                        self.blast_version = row.replace("# ", "")
                    elif "# Query: " in row:
                        self.query_header = row.replace("# Query: ", "")
                    elif "# Database: " in row:
                        self.blast_database = row.replace("# Database: ", "")
                    elif "# Fields: " in row:
                        self.fields = row.replace("# Fields: ", "").split(",")
                        self.fields = [s.strip().replace(" ", "_").replace(".", "").replace("%", "perc") for s in self.fields]
                    elif "hits found" in row:
                        self.n_hits = int(row.split(" ")[1])
                    elif "# BLAST processed" in row:
                        self.n_query = int(row.replace("# BLAST processed ", "").replace(" queries", ""))

def blastp(ref_db, query, output, word_size=6, gapopen=11, gapextend=1, evalue=0.01, outfmt=7, makedb=True):
    #if not os.path.exists(f"{ref_db}.bwt"):
    if makedb:
        subprocess.run(["makeblastdb", "-in", ref_db, "-dbtype", "prot"])
    subprocess.run(["blastp", "-query", query, "-word_size", str(word_size), "-gapopen", str(gapopen), "-gapextend", str(gapextend), "-evalue", str(evalue), "-db", ref_db, "-outfmt", str(outfmt), "-out", output])
    return os.path.realpath(output)

def blastn(ref_db, query, output, word_size=11, gapopen=5, gapextend=2, perc_identity=0.0, evalue=0.01, dust="no", outfmt=7, makedb=True):
    #if not os.path.exists(f"{ref_db}.bwt"):
    #subprocess.run(["which", "makeblastdb"])
    if makedb:
        #print(f"makeblastdb -dbtype nucl -in {ref_db}".split())
        subprocess.run(["makeblastdb", "-in", ref_db, "-dbtype", "nucl"])
        #subprocess.run(f"makeblastdb -dbtype nucl -in {ref_db}".split(), shell=True)
    #print(" ".join(["blastn", "-query", query, "-word_size", str(word_size), "-gapopen", str(gapopen), "-gapextend", str(gapextend), "-perc_identity", str(perc_identity), "-evalue", str(evalue), "-dust", dust, "-db", ref_db, "-outfmt", str(outfmt), "-out", output]))
    subprocess.run(["blastn", "-query", query, "-word_size", str(word_size), "-gapopen", str(gapopen), "-gapextend", str(gapextend), "-perc_identity", str(perc_identity), "-evalue", str(evalue), "-dust", dust, "-db", ref_db, "-outfmt", str(outfmt), "-out", output])
    return os.path.realpath(output)

def timer(func):

    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        total_time = end - start
        hour = total_time // 60 // 60
        minute = total_time // 60
        second = total_time % 60
        name = " ".join(str(func).strip("<>").split()[:2])

        if hour:
            print(f"{name} takes {hour} hour {minute} minute {second} seconds")
            return result
        elif minute:
            print(f"{name} takes {minute} minute {second} seconds")
            return result
        else:
            print(f"{name} takes {second} seconds")
            return result

    return wrapper

def check(func):

    def wrapper(*args, **kwargs):
        name = " ".join(str(func).strip("<>").split()[:2])
        print(f"TRACE: Function {name} takes in input {args} args and {kwargs} kwargs")
        result = func(*args, **kwargs)
        #print("TRACE: Function return {}".format(result))
        return result

    return wrapper
