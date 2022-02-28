#! /usr/bin/env python3

import argparse
from Bio import SeqIO
from pathlib import Path

#try:
#    from GenomicPackage.database import Database
#except ImportError:
#    import sys
#    sys.path.insert(0, "/projet/isp/galaxy_data/tools")
#    from GenomicPackage.database import Database

#def convert_accession_to_strain(database, records):
#    accession_to_strain = dict(Database.convert_accession_to_strain(database))
#    for record in records:
#        print(record.id)
#        record.id = accession_to_strain["_".join(str(record.id).split("_")[:3])]
#        yield record

#def main(fasta, database, output, translate, *args, **kwargs):
def main(fasta, format, *args, **kwargs):
    output = Path(fasta).stem + f".{format}"
    
    if format == "nexus":
        counts = SeqIO.convert(fasta, "fasta", output, format, molecule_type="DNA")
    else:
        records = SeqIO.parse(fasta, "fasta")
        #if translate:
        #    count = SeqIO.write(convert_accession_to_strain(database, records), output, "phylip")
        #    print(f"Converted {count} records")
        #else:
        count = SeqIO.write(records, output, format)
    
    print(f"Converted {count} records")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--fasta", type=str, required=True, help="Multi fasta file to convert")
    #parser.add_argument("-o", "--output", type=str, required=True, help="Output file name")
    parser.add_argument("-f", "--format", type=str, choices=["phylip", "nexus"], help="Output format")
    #parser.add_argument("-d", "--database", type=str, help="Database")
    #parser.add_argument("-t", "--translate", action="store_true", help="Translate record.id in fasta file ?")
    args = parser.parse_args()

    main(**vars(args))
