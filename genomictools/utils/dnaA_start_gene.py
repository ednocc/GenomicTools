from GenomicPackage.src import genbank
from GenomicPackage.src import search
from GenomicPackage.src import genbank_metadata as metadata

import glob
import argparse

def dnaA_start(**kwargs):
    list_gbk = glob.glob("../../genbank/*.gbk")

    i = 0
    for gbk in list_gbk:
        fna = gbk.replace('.gbk', '.fna.new')
        with open(f"{fna}", "w") as fna_out:
            GBIndex = genbank.build_index(gbk)
            seq = GBIndex.records[0].seq
            for locus in search.search_qualifier(GBIndex, genbank.q.product):
                for k2, v2 in GBIndex.index[locus].items():
                    if k2 == "CDS" and  "DnaA" in v2.data.qualifiers.get("product")[0]:
                        dnaA = v2
                        dnaA_start = dnaA.data.location.start
                        print(f"Genbank file : {gbk}\nGenome file : {fna}\ndnaA gene start at : {dnaA_start}")
                        i += 1
                        if dnaA_start != 0:
                            print("Modification ...")
                            new_seq = seq[dnaA_start:] + seq[:dnaA_start]
                            strain = metadata.extract_strain(GBIndex)
                            print(f">{strain}\n{new_seq}", file=fna_out)
                        else:
                            strain = metadata.extract_strain(GBIndex)
                            print(f">{strain}\n{seq}", file=fna_out)

        print(f"Total parsed gbk : {i}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    #parser.add_argument("", "", type=str, help="")
    args = parser.parse_args()

    dnaA_start(**vars(args))
