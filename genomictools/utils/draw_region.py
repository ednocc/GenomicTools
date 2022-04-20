#! /usr/bin/env python3

import argparse

from genomictools.genbank.Genome import Genome
from genomictools.genbank.Draw import Draw


def parse_argument():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--gbk", required=True, help="Genbank file to rotate.")
    parser.add_argument("-o", "--output", required=True, help="Output genbank file.")
    parser.add_argument("-r", "--reverse", action="store_true", help="Reverse complement the sequence.")
    parser.add_argument("-s", "--region_size", type=int, help="Reverse complement the sequence.")
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
        feature = genome.search_gene(args.gene)
        position = (feature.start + feature.end) // 2
    elif args.locustag:
        feature = genome[args.locustag]
        position = (feature.start + feature.end) // 2
    else:
        position = args.position

    region = genome.extract_region(position - args.region_size, position + args.region_size)
    print(region.features())

    title = f"Region_at_{position}"
    drawing = Draw(title, region)
    drawing.draw(args.output, "png")
