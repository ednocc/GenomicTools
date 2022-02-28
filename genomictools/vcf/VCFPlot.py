import os

import matplotlib.pyplot as plt

from genomictools.wrapper import timer, pprint
from genomitools.vcf.VCFParser import *

__all__ = ["VCFPlot"]

class VCFPlot:
    @timer
    def __init__(self, vcf):
        self.parser = VCFParser(vcf)
    
    def get_axis_data(self):
        for row in self.parser.parse_snp():
            depth = row.depth_recalibration()
            if depth: 
                self.x.append(row.qual)
                self.y.append(depth)

    def plot_row_depth_on_quality(self, output, qual_filter=20, depth_filter=5):
        self.x = [] # QUAL
        self.y = [] # Mean DP
        self.get_axis_data()
        qual_line_x = [qual_filter, qual_filter]
        qual_line_y = [0, max(self.y)]
        depth_line_x = [0, max(self.x)]
        depth_line_y = [depth_filter, depth_filter]
        plt.scatter(self.x, self.y, s=0.5, c="black")
        plt.plot(qual_line_x, qual_line_y, linewidth=0.6, color="r")
        plt.plot(depth_line_x, depth_line_y, linewidth=0.6, color="r")
        f = os.path.basename(self.parser.vcf)
        plt.title(f"SNP depth on SNP quality\n{f}")
        plt.savefig(output, dpi=800)
        plt.clf()

    #def plot_sample_depth_on_quality(self, output, qual_filter=20, depth_filter=10):
    #    qual_line_x = [qual_filter, qual_filter]
    #    qual_line_y = [0, max(self.y)]
    #    depth_line_x = [0, max(self.x)]
    #    depth_line_y = [depth_filter, depth_filter]
    #    plt.scatter(self.x, self.y, s=0.5, c="black")
    #    plt.plot(qual_line_x, qual_line_y, linewidth=0.6, color="r")
    #    plt.plot(depth_line_x, depth_line_y, linewidth=0.6, color="r")
    #    f = os.path.basename(self.parser.vcf)
    #    plt.title(f"SNP depth on SNP quality\n{f}")
    #    plt.savefig(output, dpi=800)
    #    plt.clf()

    # From roary_plots.py pangenome_matrix
    def plot_snp_along_reference(self, tree, ret=False):
        import matplotlib
        matplotlib.use('Agg')

        import matplotlib.pyplot as plt
        import seaborn as sns

        sns.set_style('white')

        import os
        import pandas as pd
        import numpy as np
        from Bio import Phylo

        t = Phylo.read(tree, 'newick')
        #t.root_with_outgroup({'name': 'PICSARLN20'})

        # Max distance to create better plots
        mdist = max([t.distance(t.root, x) for x in t.get_terminals()])

        # Load vcf
        table = self.parser.read_vcf()
        # Set index (group name)
        table.set_index('POS', inplace=True)
        table["REF"].replace('.', 0, regex=True, inplace=True)
        table.rename(columns={"REF": "GCF_000007865.1_ASM786v1"}, inplace=True)
        # Drop the other info columns - Intervalle inclusif pour columns
        table.drop(list(table.columns[:2]), axis=1, inplace=True)
        table.drop(list(table.columns[1:6]), axis=1, inplace=True)

        # Transform it in a presence/absence matrix (1/0)
        table.replace('^0:', 0, regex=True, inplace=True)
        table.replace('^1:', 1, regex=True, inplace=True)
        table.replace('^./.:', 0, regex=True, inplace=True)
        #table.replace('^./.:', np.nan, regex=True, inplace=True)

        # Sort the matrix according to tip labels in the tree
        table_sorted = table[[x.name for x in t.get_terminals()]]
        #partial_table_sorted = table[[x.name for x in t.get_terminals()]]

        #genome_size = self.parser.contig["length"]
        ## nb sample + ref
        #nb_samples = len(self.parser.samples) + 1

        #data_null = np.zeros((genome_size, nb_samples))
        #table_sorted = pd.DataFrame(data=data_null, index=range(1, genome_size + 1), columns=partial_table_sorted.columns)

        ## Put data from partial_table_sorted in table_sorted
        #table_sorted.update(partial_table_sorted)

        #if ret:
        #    return (table_sorted, partial_table_sorted)

        # Plot presence/absence matrix against the tree
        with sns.axes_style('whitegrid'):
            fig = plt.figure(figsize=(17, 10))
            #fig = plt.figure(figsize=(50, 25))

            ax1=plt.subplot2grid((1,40), (0, 10), colspan=30)
            a=ax1.matshow(table_sorted.T, cmap=plt.cm.Blues,
                            vmin=0, vmax=1,
                            aspect='auto',
                            interpolation='none',
                        )
            ax1.set_yticks([])
            ax1.set_xticks([])
            ax1.axis('off')

            ax = fig.add_subplot(1,2,1)
            # matplotlib v1/2 workaround
            try:
                ax=plt.subplot2grid((1,40), (0, 0), colspan=10, facecolor='white')
            except AttributeError:
                ax=plt.subplot2grid((1,40), (0, 0), colspan=10, axisbg='white')

            fig.subplots_adjust(wspace=0, hspace=0)

            ax1.set_title('SNP on K10')

            fsize = 12 - 0.1*table_sorted.shape[1]
            if fsize < 7:
                fsize = 7
            with plt.rc_context({'font.size': fsize}):
                Phylo.draw(t, axes=ax, 
                        show_confidence=False,
                        label_func=lambda x: str(x)[:10],
                        xticks=([],), yticks=([],),
                        ylabel=('',), xlabel=('',),
                        xlim=(-mdist*0.1,mdist+mdist*0.45-mdist*table_sorted.shape[1]*0.001),
                        axis=('off',),
                        title=('Tree\n(%d strains)'%table_sorted.shape[1],), 
                        do_show=False,
                        )

            plt.savefig('SNP matrix.png', dpi=600)
            plt.clf()
