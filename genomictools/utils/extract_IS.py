#! /usr/bin/env python3

import argparse
import sys
from Bio import SeqIO
#from Bio.SeqRecord import SeqRecord

try:
    from GenomicPackage.genbank import Genome
except ImportError:
    sys.path.insert(0, "/projet/isp/galaxy_data/tools")
    from GenomicPackage.genbank import Genome

def main(args):
    GBIndex = Genome(args.genbank)

    print("")
    print(GBIndex.path)
    print("")
    
    element = {}

    i = 0
    for locustag, feature in GBIndex:
        #if not GBIndex.index[l]['CDS'].data.qualifiers.get('pseudo', []):
        # term = family transposase if you search for IS element annotation
        if args.term in feature.product:
            r = [feature.product, str(feature.start), str(feature.end)]
            # Add 1 if IS is already in the dictonary keys
            elem = r[0].split(" ")[0]
            if elem in element.keys():
                element[elem] += 1
            else:
                element[elem] = 1

            print("\t".join(r))
            #seq = feature.seq
            record = feature.to_seqrecord()
            SeqIO.write(record, args.output, "fasta")
            #print(f">{r[1]}_{r[2]}\n{seq}", file=args.output)
            i += 1

    print("\n#### stat ####\n")

    for k, v in element.items():
        print(f"{k} : {v}")

    print(f"Total : {i}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Extract position of IS110 family insertion sequence")
    parser.add_argument("-g", "--genbank", required=True, type=str, help="genbank file of the strain")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout, help="Output file (default = stdout)")
    parser.add_argument("-t", "--term", required=True, type=str, help="search term")

    args = parser.parse_args()

    main(args)
