import sys
import os
import itertools
import shutil
import subprocess
from datetime import datetime

from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.units import cm

from Bio import SearchIO

from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import CrossLink

__all__ = ["Draw"]

class Draw:
    def __init__(self, title, *genome_region):
        self.tmp_dir = Path(__file__).parent / f"tmp_{datetime.now().strftime('%d%m%y_%H%M%S_%f')}"
        if not self.tmp_dir.exists():
            self.tmp_dir.mkdir()

        if len(genome_region) == 0:
            raise ValueError(f"Nothing to draw.")
        elif len(genome_region) > 1:
            # Upper track has the biggest number
            self.genome_region = list(reversed(genome_region))
        else: # genome_region == 1
            self.genome_region = genome_region

        self.diagram = GenomeDiagram.Diagram(title)
        self.track = {}

        for n, region in enumerate(self.genome_region, 1):
            track = self.diagram.new_track(n, name=region.gbk, scale=0)
            track_feature_set = track.new_set()
            for locustag, feature in region:
                print(feature.locustag, feature.gene, feature.start, feature.end, sep="\t")
                feat_sigil = "ARROW"
                feat_arrowshaft = 0.6
                feat_border_color = colors.black
                feat_color = colors.white
                track_feature_set.add_feature(feature.feature, sigil=feat_sigil, color=feat_color, border=feat_border_color, arrowshaft_height=feat_arrowshaft)
            self.track[region.records[0].id] = [region, track, track_feature_set]


    def __del__(self):
        shutil.rmtree(self.tmp_dir)


    def _pairwise(self, iterable):
        # itertools.pairwise substitute as this function is introduces in python version 3.10
        # pairwise('ABCDEFG') --> AB BC CD DE EF FG
        a, b = itertools.tee(iterable)
        next(b, None)
        return zip(a, b)


    def blast(self):
        # itertools.pairwise('ABCDEFG') --> AB BC CD DE EF FG
        for (region1_id, (region1, track1, feature_set1)), (region2_id, (region2, track2, feature_set2)) in self._pairwise(self.track.items()):
            region1_fasta = f"{self.tmp_dir}/{region1.gbk.replace('.gbk', '_' + str(region1.left_origin) + '_' + str(region1.right_origin) + '.fasta')}"
            region1.format(region1_fasta, "fasta")
            region2_fasta = f"{self.tmp_dir}/{region2.gbk.replace('.gbk', '_' + str(region2.left_origin) + '_' + str(region2.right_origin) + '.fasta')}"
            region2.format(region2_fasta, "fasta")

            query = region1_fasta
            database = region2_fasta
            blast_result = f"{self.tmp_dir}/{track1.name}_on_{track2.name}.blast"

            # If index database file doesn't exists or database file if newer than blast index file then rebuild blast database
            if not os.path.exists(database + ".nhr") or os.path.getmtime(database) > os.path.getmtime(database + ".nhr"):
                print(f"Creating BLAST database from input ({database})", file=sys.stderr)
                subprocess.run(f"makeblastdb -in {database} -dbtype nucl 1>/dev/null", shell=True)
            else:
                print(f"BLAST database for {database} already exists", file=sys.stderr)
            print(f"Performing BLAST of ({query}) against input {database}", file=sys.stderr)
            subprocess.run(f"blastn -query {query} -db {database} -outfmt 7 -out {blast_result}", shell=True)

            print(f"Parsing BLAST results from {blast_result}", file=sys.stderr)

            # Parse blast results
            # QueryResult   ->    [Hits]    ->      [Hsps]
            # per query         per results         per alignment
            #                   in database/
            #                   reference
            # Only 1 hit with many alignment (hsp)
            for query in SearchIO.parse(blast_result, "blast-tab", comments=True):
                #query = record.id
                #primer_name, primer_id = query.split("_")
                for i, hit in enumerate(query.hits, 1):
                    for n, hsp in enumerate(hit.hsps, 1):
                        query_link_coordinate = (self.track[query.id][1], *hsp.query_range)
                        hit_link_coordinate = (self.track[hit.id][1], *hsp.hit_range)

                        link = CrossLink(query_link_coordinate, hit_link_coordinate, colors.firebrick, colors.lightgrey)
                        self.diagram.cross_track_links.append(link)


    def draw(self, prefix, output_format):
        if len(self.genome_region) > 1:
            self.blast()

        self.diagram.draw(format="linear", pagesize=(1920, 250 * ((len(self.genome_region) * 2) - 1)), fragments=1, fragment_size=0.99, track_size=0.2, x=0.01, y=0.01)#, start=0, end=15903, x=0.01, y=0.01)
        self.diagram.write(prefix + "." + output_format, output_format)


if __name__ == '__main__':
    from Genome import Genome

    h37Rv_genome = Genome(sys.argv[1])
    g20x_contigs = Genome(sys.argv[2])

    katG_h37rv = h37Rv_genome.extract_region(2151422, 2162770)
    katG_g20x = g20x_contigs.extract_region(1110, 8800, reverse=True)

    drawing = Draw("katG region", katG_h37rv, katG_g20x)
    drawing.draw("katG_region", "png")
