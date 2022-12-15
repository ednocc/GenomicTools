#! /usr/bin/env python3

import argparse
from genomictools.genbank.Genome import Genome


def parse_argument():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--gbk", required=True, help="Genbank file to rotate.")
    parser.add_argument("-o", "--output", required=True, help="Output genbank file.")
    parser.add_argument("-r", "--reverse", action="store_true", help="Reverse complement the sequence.")
    parser.add_argument("--start", type=int, help="Position in 0-based coordinate which will serve as the beginning of the region to extract.")
    parser.add_argument("--stop", type=int, help="Position in 0-based coordinate which will serve as the end of the region to extract. Note that this position is exclude from the region")
    parser.add_argument("--record_id", type=str, default="",  help="Chromosome or contig accession (name in the LOCUS line). [default = first chromosome or contig]")
    #group = parser.add_mutually_exclusive_group(required=True)
    #group.add_argument("-g", "--gene", help="Gene name which will serve as the new origin.")
    #group.add_argument("-l", "--locustag", help="Locustag which will serve as the new origin.")
    #group.add_argument("-p", "--position", type=int, help="Position in 0-based coordinate which will serve as the new origin.")

    parsed_args = parser.parse_args()

    if not parsed_args.record_id:
        parsed_args.record_id = 0

    return parsed_args


def main():
    args = parse_argument()

    genome = Genome(args.gbk)
    new_genome = genome.extract_region(args.start, args.stop, record_id=args.record_id, reverse=args.reverse)

    new_genome.format(args.output, "genbank")


if __name__ == '__main__':
    main()
