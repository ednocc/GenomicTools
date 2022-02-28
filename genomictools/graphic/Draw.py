#!/usr/bin/env python
# coding: utf-8

# In[89]:


import matplotlib.pyplot as plt
from matplotlib.path import Path as _path
import matplotlib.patches as _patches


# In[ ]:


#verts = [
#    (0., 0.), # left, bottom
#    (0., 1.), # left, top
#    (1., 1.), # right, top
#    (1.5, 0.5),
#    (1., 0.), # right, bottom
#    (0., 0.), # ignored
#    ]
#
#codes = [Path.MOVETO,
#         Path.LINETO,
#         Path.LINETO,
#         Path.LINETO,
#         Path.LINETO,
#         Path.CLOSEPOLY,
#         ]
#
#path = Path(verts, codes)
#
#fig = plt.figure()
#ax = fig.add_subplot(111)
#ax.axis([-1, 2010, -1, 2])
#patch = patches.PathPatch(path, facecolor='red', lw=0)
#ax.add_patch(patch)
##ax.set_xlim(-2,2)
##ax.set_ylim(-2,2)


# In[199]:


class Feature:
    height = 1
    init = 0.
    
    def __init__(self, start=0, stop=1, strand=None, color="auto", label="", **kwargs):
        self.start = float(start)
        self.stop = float(stop)
        self.strand = strand
        self.draw_options = kwargs
        
        if self.strand:
            if self.strand == -1 and start > stop or self.strand == 1 and start > stop:
                self.start, self.stop = self.stop, self.start
        else:
            if start > stop:
                self.start, self.stop = self.stop, self.start
                
        self.len = self.stop - self.start
        
        if color == "auto":
            if self.strand == 1:
                self.color = "pink"
            elif strand == -1:
                self.color = "lightblue"
            else:
                self.color = "grey"
        else:
            self.color = color

        self.label = label


class Arrow(Feature):
    arrow = 25.
    
    codes = [_path.MOVETO,
         _path.LINETO,
         _path.LINETO,
         _path.LINETO,
         _path.LINETO,
         _path.CLOSEPOLY,
         ]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        if self.strand == 1:
            self.verts = [
                (self.start, Feature.init), # left, bottom
                (self.start, Feature.height), # left, top
                (self.start + self.len - Arrow.arrow, Feature.height), # right, top
                (self.start + self.len, Feature.height / 2), # arrow
                (self.start + self.len - Arrow.arrow, Feature.init), # right, bottom
                (self.start, Feature.init), # ignored return left, bottom
                ]
        elif self.strand == -1:
            self.verts = [
                (self.stop, Feature.init), # right, top
                (self.stop, -Feature.height), # right, bottom
                (self.stop - self.len + Arrow.arrow, -Feature.height), # left, bottom
                (self.stop - self.len, -Feature.height / 2), # arrow
                (self.stop - self.len + Arrow.arrow, Feature.init), # left, top
                (self.stop, Feature.init), # ignored return right, top
                ]
        else: # strand = None
            self.verts = [
                (self.start, -Feature.height), # left, bottom
                (self.start, Feature.height), # left, top
                (self.start + self.len - Arrow.arrow, Feature.height), # right, top
                (self.start + self.len, Feature.init), # arrow
                (self.start + self.len - Arrow.arrow, -Feature.height), # right, bottom
                (self.start, -Feature.height), # ignored return left, bottom
                ]
        
        self.path = _path(self.verts, self.__class__.codes)
        self.patch = _patches.PathPatch(self.path, facecolor=self.color, alpha=0.5, lw=5)


class Box(Feature):
    codes = [_path.MOVETO,
         _path.LINETO,
         _path.LINETO,
         _path.LINETO,
         _path.CLOSEPOLY,
         ]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        if self.strand == 1:
            self.verts = [
                (self.start, Feature.init), # left, bottom
                (self.start, Feature.height), # left, top
                (self.start + self.len, Feature.height), # right, top
                (self.start + self.len, Feature.init), # right, bottom
                (self.start, Feature.init), # ignored return left, bottom
                ]
        elif self.strand == -1:
            self.verts = [
                (self.stop, Feature.init), # right, top
                (self.stop, -Feature.height), # right, bottom
                (self.stop - self.len, -Feature.height), # left, bottom
                (self.stop - self.len, Feature.init), # left, top
                (self.stop, Feature.init), # ignored return right, top
                ]
        else: # strand = None
            self.verts = [
                (self.start, -Feature.height), # left, bottom
                (self.start, Feature.height), # left, top
                (self.start + self.len, Feature.height), # right, top
                (self.start + self.len, -Feature.height), # right, bottom
                (self.start, -Feature.height), # ignored return left, bottom
                ]
        
        self.path = _path(self.verts, self.__class__.codes)
        self.patch = _patches.PathPatch(self.path, facecolor=self.color, alpha=0.3, lw=5)


class Graph:
    def __init__(self, feats, dpi=300):
        self.dpi = dpi
        self.feats = feats
        self.min_x = min([f.start for f in self.feats]) - 1
        self.max_x = max([f.stop for f in self.feats]) + 1
        
        self.fig = plt.figure(figsize=(50, 20))
        self.ax = self.fig.add_subplot(111)
        self.ax.axis([self.min_x, self.max_x, -2, 3])
        
        self.ax.spines["top"].set_color("none")
        #self.ax.spines["bottom"].set_color("none")
        self.ax.spines["right"].set_color("none")
        self.ax.spines["left"].set_color("none")
        self.ax.yaxis.set_ticks([])
        self.ax.yaxis.set_ticklabels([])
        #self.ax.xaxis.set_ticks([])
        #self.ax.xaxis.set_ticklabels([])
        
        self.fontsize = 36
        
    def draw(self, output=None):
        for feat in self.feats:
            self.ax.add_patch(feat.patch)
            #print(feat.label, feat.strand)

            if not feat.draw_options.get("rotation", None):
                feat.draw_options["rotation"] = 45 if feat.strand == 1 or not feat.strand else -45

            if not feat.draw_options.get("fontsize", None):
                feat.draw_options["fontsize"] = self.fontsize

            if not feat.draw_options.get("weight", None):
                feat.draw_options["weight"] = "bold"

            if feat.strand == 1 or not feat.strand:
                xy = (feat.stop - (feat.len / 2), 1.25)
                self.ax.annotate(feat.label, xy=xy, **feat.draw_options)
            else:
                xy = (feat.stop - (feat.len / 2), -1.5)
                self.ax.annotate(feat.label, xy=xy, **feat.draw_options)

        if output:
            plt.savefig(output, dpi=self.dpi)
        else:
            plt.show()


if __name__ == '__main__':
    feat1 = Box(start=3867288 - 2000, stop=3867288, label="A")
    feat2 = Box(start=3867288, stop=3868741, label="15")
    feat3 = Arrow(start=3867288 + 226, stop=3868741 - 14, strand=1, label="transposase", fontsize=10)
    feat4 = Box(start=3868741, stop=3868741 + 2000, label="B")
    feats = [feat1, feat2, feat3, feat4]

    #feat1 = Arrow(start=0, stop=1000, label="LOCUS_1")
    ##feat2 = Arrow(start=500, stop=1500, strand=-1, label="LOCUS_2")
    ##feat3 = Arrow(start=2500, stop=2000, strand=1, label="LOCUS_3")
    ##feat4 = Arrow(start=4000, stop=500, strand=-1, label="LOCUS_4")
    #feat5 = Arrow(start=3250, stop=2500, strand=-1, label="LOCUS_5")
    #feat6 = Box(start=1600, stop=1800, strand=-1, label="LOCUS_6")
    #feat7 = Box(start=2000, stop=2200, label="LOCUS_7")
    #feats = [feat1, feat6, feat7, feat5]
    #feats = [feat6, feat7, feat2, feat1, feat3, feat4, feat5]
    
    #feats_pos = [Feature(i, i + 1000, label=f"LOCUS_{i}") for i in range(0, 100000, 1500)]
    #feats_neg = [Feature(i, i + 1000, strand=-1, label=f"LOCUS{i}") for i in range(0, 100000, 1500)]
    #feats = feats_pos + feats_neg
    
    graph = Graph(feats)
    graph.draw("features.png")

