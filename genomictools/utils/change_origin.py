import argparse
from genomictools.genbank.Genome import Genome


def parse_argument():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--gbk", required=True, help="Genbank file to rotate.")
    parser.add_argument("-o", "--output", required=True, help="Output genbank file.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-l", "--locustag", help="Locustag which will serve as the new origin.")
    group.add_argument("-g", "--gene", help="Gene name which will serve as the new origin.")

    parsed_args = parser.parse_args()
    return parsed_args

if __name__ == '__main__':
    args = parse_argument()

    genome = Genome(args.gbk)
    if args.gene:
        new_genome = genome.change_origin(gene_name=args.gene)
    else:
        new_genome = genome.change_origin(locustag=args.locustag)
    new_genome.format(args.output, "genbank")
