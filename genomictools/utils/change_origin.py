#! /usr/bin/env python3

import argparse
from genomictools.genbank.Genome import Genome


def parse_argument():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--gbk", required=True, help="Genbank file to rotate.")
    parser.add_argument("-o", "--output", required=True, help="Output genbank file.")
    parser.add_argument("-r", "--reverse", action="store_true", help="Reverse complement the sequence.")
    parser.add_argument("-d", "--debug", action="store_true", help="Enable debugging messages.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-g", "--gene", help="Gene name which will serve as the new origin.")
    group.add_argument("-l", "--locustag", help="Locustag which will serve as the new origin.")
    group.add_argument("-p", "--position", type=int, help="Position which will serve as the new origin.")

    parsed_args = parser.parse_args()
    return parsed_args

if __name__ == '__main__':
    args = parse_argument()

    genome = Genome(args.gbk)
    if args.gene:
        new_genome = genome.change_origin(gene_name=args.gene, reverse=args.reverse, debug=args.debug)
    elif args.locustag:
        new_genome = genome.change_origin(locustag=args.locustag, reverse=args.reverse, debug=args.debug)
    else:
        new_genome = genome.change_origin(position=args.position, reverse=args.reverse, debug=args.debug)

    new_genome.format(args.output, "genbank")
