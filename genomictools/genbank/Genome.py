import os
import re

from pathlib import Path
from collections import OrderedDict

from multiprocessing import Pool

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC

from reportlab.lib import colors
from reportlab.lib.units import cm

from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import CrossLink

#pip install dna_features_viewer
from dna_features_viewer import GraphicFeature, GraphicRecord
import matplotlib.pyplot as plt

from genomictools.wrapper import pprint, timer

#import pickle
#from GenomicPackage.hashfile import HashFile

from genomictools.genbank.Feature import Feature

#try:
#    from Feature import Feature
#except:
#    from .Feature import Feature
#from .Metadata import Metadata

import itertools
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor, as_completed

FORWARD_STRAND = 1
REVERSE_STRAND = -1

__all__ = ["FORWARD_STRAND", "REVERSE_STRAND", "Genome", "PartialGenome", "process_multiple_genome", "parallel_genomes", "sequential_genomes"]


class BaseGenome:
    def __init__(self, genbank):
        self.gbk = os.path.basename(genbank)
        self.path = os.path.realpath(genbank)
        self.fasta = None
        self.extension = os.path.splitext(self.path)[1]
        self.name = "_".join(self.gbk.strip(self.extension).split("_")[:3])
        self.major_type = ["CDS", "tRNA", "rRNA", "ncRNA", "tmRNA"]
        self.records = []
        self.seq = []
        self.qualifiers = OrderedDict()
        self.other_features = []
        self.index = OrderedDict()
        self.keys_index = []
        self.strain_sep = "-"

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return f"""
File    : {self.path}
"""

    def __iter__(self):
        return GenomeIterator(self)

    def __len__(self):
        length = 0
        for record in self.index.keys():
            length += len(record.keys())
        return length

    def __eq__(self, other):
        state = True
        if self.index != other.index:
            state = False
        return state

    def __bool__(self):
        return True if self.index else False

    def __getitem__(self, key):
        if not self:
            return None

        if isinstance(key, slice):
            try:
                if self[key.start].contig_id != self[key.stop].contig_id:
                    return None
                else:
                    contig = self[key.start].contig_id

                if self.keys_index.index(key.start) > self.keys_index.index(key.stop):
                    start = key.stop
                    stop = key.start
                else:
                    start = key.start
                    stop = key.stop
            except ValueError as e:
                print(f"KeyError : {e}")


            between = OrderedDict()
            save = False
            for locustag, feature in self.index.items():
                if locustag == start:
                    save = True
                    left_origin = feature.start
                if save:
                    if locustag == stop:
                        save = False
                        right_origin = feature.end
                    between[locustag] = feature

            seqrecord_features = [self.records[contig].features[0]] # Source
            seqrecord_features += [v.lshift(left_origin).feature for k, v in between.items()] # SeqFeature
            seqrecord = [SeqRecord(self.seq[contig][int(left_origin):int(right_origin)], 
                                        id=self.records[contig].id, 
                                        name=self.records[contig].name, 
                                        description=self.records[contig].description, 
                                        dbxrefs=self.records[contig].dbxrefs, 
                                        features=seqrecord_features, 
                                        annotations=self.records[contig].annotations, 
                                        letter_annotations=self.records[contig].letter_annotations)]

            return PartialGenome(self.gbk, between, seqrecord, left_origin, right_origin)
        else:
            if isinstance(key, int):
                try:
                    key = self.keys_index[key]
                except:
                    return None

            return self.index.get(key, None)
            
    def parse_records(self):
        for record in SeqIO.parse(self.path, "genbank"):
            self.records.append(record)

    def parse_seq(self):
        for record in self.records:
            self.seq.append(record.seq)

    #def get_metadata(self):
    #    self.metadata = Metadata(self)

    #def update_keys_index(self):
    #    self.keys_index = list(self.index.keys())
        
    def show_all_qualifiers(self):
        pprint(self.qualifiers)

    def format(self, output, fmt):
        count = SeqIO.write(self.records, output, fmt)
        return output

    def to_fasta(self, output):
        return self.format(output, "fasta")

    def check_fasta(self, output=""):
        if not output:
            self.fasta = self.path.replace(self.extension, '.fna')
        else:
            self.fasta = output

        if not os.path.exists(self.fasta):
            self.to_fasta(self.fasta)
        
        return self.fasta

    #####  Metadata #####

    @property
    def species(self):
        species = self.records[0].features[0].qualifiers.get("organism")[0]
        if "subsp." in species:
            return " ".join(species.split()[:4])
        else:
            return " ".join(species.split()[:2])
    
    @property
    def strain(self):
        name_fields = ["strain", "isolate"]
        for field in name_fields:
            strain_name = self.records[0].features[0].qualifiers.get(field, None)
            if strain_name:
                return strain_name[0].replace(" ", self.strain_sep).replace("/", self.strain_sep)
        
        raise ValueError("No field name found.")
        #return self.source.qualifiers.get("strain")[0].replace(" ", self.strain_sep).replace("/", self.strain_sep)
        
    @property
    def short_strain(self):
        #return self.source.qualifiers.get("strain")[0].replace(" ", "").replace("/", "").replace("-", "").replace("_", "")
        return self.strain.replace("-", "").replace("_", "")
        
    @property
    def dbxref(self):
        return self.records[0].features[0].qualifiers.get("db_xref", "NA")[0]

    @property
    def species_code(self):
            return self.species.split()[-1][:4]

    def extract_sequence_id(self):
        return [record.id for record in self.records] 
    
    def extract_plasmid(self):
        if len(self.records) > 1:
            return [p.id for p in self.records[1:]]
        elif len(self.records) > 20:
            raise ValueError("Your genome seems to be an assembly. This method can be use only on complete genome.")
        else:
            return []
    
    def count_gene(self):
        i = [0 for n in range(0, len(self.records))]
        for n in range(0, len(self.records)):
            for locustag, feature in self:
                if feature.contig_id == n:
                    i[n] += 1
        return i
    
    def count_tRNA(self):
        i = [0 for n in range(0, len(self.records))]
        for n in range(0, len(self.records)):
            for locustag, feature in self:
                if feature.contig_id == n and feature.type == 'tRNA':
                    i[n] += 1
        return i
    
    def count_rRNA(self):
        i = [0 for n in range(0, len(self.records))]
        for n in range(0, len(self.records)):
            for locustag, feature in self:
                if feature.contig_id == n and feature.type == 'rRNA':
                    i[n] += 1
        return i
    
    def count_coding_and_pseudo(self):
        i = [[0, 0] for n in range(0, len(self.records))]
        for n in range(0, len(self.records)):
            for locustag, feature in self:
                if feature.contig_id == n and feature.type == 'CDS' and not feature.pseudo:
                    i[n][0] += 1
                elif feature.contig_id == n and feature.pseudo:
                    i[n][1] += 1
                else:
                    continue
        return i
    
    def count_hypothetical(self):
        i = [0 for n in range(0, len(self.records))]
        for n in range(0, len(self.records)):
            for locustag, feature in self:
                if feature.contig_id == n and feature.product == "hypothetical protein":
                    i[n] += 1
        return i
    
    def calc_GC_content(self):
        # 6 --> arbitrary cutoff (plasmid in genbank of complete genome)
        # Draft genome currently have more than 6 contigs
        if len(self.records) < 7:
            return round(GC(self.seq[0]), 2)
        else:
            pGC = 0
            i = 0
            for s in self.seq:
                pGC += GC(s)
                i += 1
            pGC = pGC // i
            return round(pGC, 2)


class GenomeIterator:
    def __init__(self, genome):
        self.genome = genome
        self.contig = 0
        self.n_features_contig = 0
        self.current_feature = 1
        self.cursor = 0
        self.max = len(genome) # Total number of locustag in all records

    def __iter__(self):
        return self

    def __next__(self):
        if self.cursor >= self.max: # End of the file - StopIteration
            raise StopIteration
        self.n_features_contig = len(genome.index[self.contig])
        locustag = list(self.genome.index[self.contig].keys())[self.current_feature]
        for feat_id in self.genome.index[self.contig][locustag]:
            feature = genome.records[self.contig].features[feat_id]
            if feature.type in genome.major_type: # Feature found - Next locustag
                self.cursor += 1
                self.current_feature += 1
                if self.current_feature >= self.n_features_contig: # End of this record - Next record
                    self.contig += 1
                    self.current_feature = 0
                return (locustag, feature)


class Genome(BaseGenome):
    #@timer
    def __init__(self, genbank):
        super().__init__(genbank)
        #self.cached = False
        #self.md5hasher = HashFile()
        #self.check_cached_data()
        self._dict_index()
        #self.update_keys_index()
        #self.source = self.records[0].features[0]

    def __repr__(self):
        return f"""
Strain  : {self.strain}
Species : {self.species}
File    : {self.path}
GC%     : {self.calc_GC_content()}
"""

    def _dict_index(self):
        for (contig_id, record) in enumerate(SeqIO.parse(self.path, "genbank")):
            self.records.append(record)
            self.index[contig_id] = {}
            for (index, feature) in enumerate(record.features[1:]):
                try:
                    locustag = feature.qualifiers.get("locus_tag")[0]
                except TypeError as e:
                    continue
                if not locustag in self.index[contig]:
                    self.index[contig_id][locustag] = [index]
                else:
                    self.index[contig_id][locustag].append(index)

    #def _parse_qualifiers(self, feature_type, qualifiers):
    #def _parse_qualifiers(self, feature):
    #    for qualifier in feature.qualifiers.items():
    #        if feature.type in self.qualifiers.keys():
    #            if qualifier[0] not in self.qualifiers[feature.type]:
    #                self.qualifiers[feature.type].append(qualifier[0])
    #        else:
    #            self.qualifiers[feature.type] = [qualifier[0]]

    ###### Search #######

    def search_gene(self, gene):
        #result = OrderedDict()
        for locustag, feature in self:
            if gene in feature.gene:
                #result[locustag] = feature
                return feature
        #return PartialGenome(self.gbk, result, self.records)
    
    def search_feature(self, type="CDS"):
        result = OrderedDict()
        if type in self.major_type:
            for locustag, feature in self:
                if feature.type == type:
                    result[locustag] = feature
        return PartialGenome(self.gbk, result, self.records)

    def search_product(self, product, regex=False):
        result = OrderedDict()
        for locustag, feature in self:
            if regex:
                s = re.search(product, feature.product)
                if s:
                    result[locustag] = feature
            else:
                if product in feature.product:
                    result[locustag] = feature
        return PartialGenome(self.gbk, result, self.records)
    

    def search_position(self, pos, contig=0):
        contig = self.translate_contigs_name(contig)

        try:
            for locustag, feature in self:
                if feature.contig_id == contig:
                    if pos in feature:
                        # If feature len is the same as contig len so this is probably an error
                        if feature.end == len(self.records[contig]):
                            continue
                        return (locustag, feature)
                    else:
                        if pos < feature.start:
                            return ("Intergenic", None)
            else:
                return ("Intergenic", None)
        except Exception as e:
            raise ValueError(f"contig params must be int or str. {type(contig)}.\n({e})")

    # TODO : manage fragmented genome
    def extract_locustag_before_after(self, start, end, contig=0, pseudo=True):
        before = []
        after = []
        #contig = self.translate_contigs_name(contig)

        try:
            for locustag, feature in self:
                if feature.start < start:
                    before.append((locustag, feature))
                if feature.start > end:
                    after.append((locustag, feature))
        except Exception as e:
            raise ValueError(f"contig params must be int or str. {type(contig)}.\n({e})")

        # Filtering
        if not pseudo:
            before = [(l, f) for l, f in before if not f.pseudo]
            after = [(l, f) for l, f in after if not f.pseudo]

        ret = OrderedDict([before[-1], after[0]])
        return PartialGenome(self.gbk, ret, self.records)

    def extract_locustag_in_region(self, left, right, contig=0):
        return self.extract_region(left, right, contig)

    def extract_dna_in_region(self, left, right, contig=0):
        contig = self.translate_contigs_name(contig)
        return self.extract_region(left, right, contig).seq[0]

    def translate_contigs_name(self, name):
        if isinstance(name, str):
            for contig, record in enumerate(self.records):
                if record.id == name or record.name == name:
                    return contig
            else:
                raise Exception(f"Contig name or id was not be found in this genome. {contig}")
        else:
            if name >= len(self.records):
                raise ValueError("Contigs value is too high. Given: {name} Max: {len(self.records)}")
            return int(name)

    def extract_region(self, left, right, contig=0, strict=True, reverse=False):
        between = OrderedDict()
        other_features = []
        contig = self.translate_contigs_name(contig)

        #if contig is not None:
        try:
            for locustag, feature in self:
                if feature.contig_id == contig:
                    #if feature.strand == 1:
                    #    start = feature.start
                    #    end = feature.end
                    #else:
                    #    # Inversion du start / end pour conserver la même condition
                    #    start = feature.end
                    #    end = feature.start
                    #if (start <= left and left <= end) or (left <= start and end <= right) or (start <= right and right <= end):
                    if feature.start <= left <= feature.end: # overlapping feature on the left border
                        if not strict:
                            left = feature.start
                        between[locustag] = feature
                    elif left <= feature.start and feature.end <= right: # feature intern to the region
                        between[locustag] = feature
                    elif feature.start < right and right <= feature.end:
                        if not strict: # overlapping feature on the right border
                            right = feature.end
                        between[locustag] = feature
                    else:
                        continue

            for feature in self.other_features:
                if feature.contig_id == contig:
                    if feature.start <= left <= feature.end: # overlapping feature on the left border
                        other_features.append(feature)
                    elif left <= feature.start and feature.end <= right: # feature intern to the region
                        other_features.append(feature)
                    elif feature.start < right and right <= feature.end:
                        other_features.append(feature)
                    else:
                        continue

            seqrecord_features = [self.records[contig].features[0]] # Source
            if reverse:
                new_seq = self.seq[contig][int(left):int(right)].reverse_complement()
                between = {locustag: feature.lshift(int(left)).reverse_complement(len(new_seq)) for locustag, feature in zip(reversed(list(between.keys())), reversed(list(between.values())))}
                other_features = [feature.lshift(int(left)).reverse_complement(len(new_seq)) for feature in reversed(other_features)]
                seqrecord_features += [feature.feature for feature in list(between.values()) + other_features] # SeqFeature
            else:
                new_seq = self.seq[contig][int(left):int(right)]
                between = {locustag: feature.lshift(int(left)) for locustag, feature in between.items()}
                other_features = [feature.lshift(int(left)) for feature in other_features]
                seqrecord_features += [feature.feature for feature in list(between.values()) + other_features] # SeqFeature

            seqrecord = [SeqRecord(new_seq, 
                                        #id=f"{self.records[contig].id}_{int(left)}:{int(right)}", 
                                        id=self.records[contig].id, 
                                        name=self.records[contig].name, 
                                        description=self.records[contig].description, 
                                        dbxrefs=self.records[contig].dbxrefs, 
                                        features=seqrecord_features, 
                                        annotations=self.records[contig].annotations, 
                                        letter_annotations=self.records[contig].letter_annotations)]
        except Exception as e:
            raise TypeError(f"contig params must be int (contig index) or str (contigs accession). {type(contig)}.\n({e})")
        

        return PartialGenome(self.gbk, between, other_features, seqrecord, left, right)

    def dnaA_at_origin(self, output=""):
        if output:
            new_fna = Path(self.check_fasta(output))
        else:
            new_fna = Path(self.check_fasta() + ".tmp")

        dnaA_locustag, dnaA_feature = [(locustag, feature) for locustag, feature in self.search_product("DnaA")]
        if not type(dnaA_feature) == list:
            dnaA_start = dnaA_feature.start
            if dnaA_start not in range(0, 11):
                self.seq[0] = self.seq[0][dnaA_start:] + self.seq[0][:dnaA_start]
                self.to_fasta(new_fna)
                return new_fna
            else:
                return self.fasta
        else:
            raise TypeError(f"'dnaA_feature' variable don't have the right type : {type(dnaA_feature)}")
        

    def change_origin(self, locustag=None, gene_name=None):
        # Find the locustag or gene_name
        if gene_name:
            starting_feature = self.search_gene(gene_name)
            start = starting_feature.start
        if locustag:
            starting_feature = self[locustag]
            start = starting_feature.start

        seg_right = self.extract_region(start, len(self.seq[0]))
        seg_left = self.extract_region(0, start)
        #Debug
        #seg_left.format("Test_left.gbk", "genbank")
        #seg_right.format("Test_right.gbk", "genbank")
        new_genome = seg_right + seg_left
        return new_genome


    #def check_cached_data(self):
    #    md5hash = self.md5hasher.hash_file(self.gbk)
    #    path_cache_dir = Path(__file__).resolve().parent.joinpath(".cache")
    #    if not path_cache_dir.is_dir():
    #        path_cache_dir.mkdir()
    #    cache_list = [f.name for f in path_cache_dir.glob("*")]
    #    path_cache_file = path_cache_dir.joinpath(md5hash)
    #    if md5hash in cache_list:
    #        self.cached = True
    #        with open(path_cache_file, "rb") as fd:
    #            self.index = pickle.load(fd)
    #    else:
    #        # Save index object in binary file named with the md5 hash of the genbank file used to created it
    #        self._dict_index()
    #        with open(path_cache_file, "wb") as fd:
    #            pickle.dump(self.index, fd)


class PartialGenome(BaseGenome):
    def __init__(self, genbank, index, other_features, seqrecord, left_origin=0, right_origin=0):
        super().__init__(genbank)
        self.left_origin = left_origin
        self.right_origin = right_origin
        self.records = seqrecord
        self.source = self.records[0].features[0]
        self.parse_seq()
        self.index = index
        self.other_features = other_features
        #self.update_keys_index()


    def __add__(self, other):
        if not isinstance(other, PartialGenome):
            raise ValueError(f"other must be instance of PartialGenome, found {type(other)}")
        if len(self.records) > 1:
            raise ValueError(f"To concatenate two PartialGenome, there must be only one sequence (number of sequence {len(self.seq)}).")

        #print(self.records[0].features[1])
        #print(other.records[0].features[1])
        #print(other.records[0].features[-1])

        # Construct new seq
        new_seq = self.seq[0] + other.seq[0]

        # Construct new index by adding length of self.seq[0] to start of each feature from other
        new_index = self.index.copy()
        new_index.update({locustag: feature + len(self.seq[0]) for locustag, feature in other.index.items()})
        new_other_features = self.other_features.copy() + other.other_features.copy()
        seqrecord_features = [self.records[0].features[0]] # Source
        seqrecord_features += [feature.feature for feature in list(new_index.values()) + new_other_features] # SeqFeature
        #print(seqrecord_features[1])

        # Construct new seqrecord based on new_index and new_seq
        seqrecord = [SeqRecord(new_seq, 
                                    id=self.records[0].id, 
                                    name=self.records[0].name, 
                                    description=self.records[0].description, 
                                    dbxrefs=self.records[0].dbxrefs, 
                                    features=seqrecord_features, 
                                    annotations=self.records[0].annotations, 
                                    letter_annotations=self.records[0].letter_annotations)]

        # Return new PartialGenome instance
        return PartialGenome(self.gbk, new_index, new_other_features, seqrecord)


    def reset_origin_coordinate(self):
        for locustag, feature in self:
            feature = feature + self.left_origin

    def draw(self, output, output_format, diagram=None, feature_set=None, diagram_name="MyDiagram", format="linear", pagesize=(1920, 80), key_color=None, key_color_args=None):
        self.reset_origin_coordinate()

        if not diagram and not feature_set:
            diagram = GenomeDiagram.Diagram(diagram_name)
            track_for_features = diagram.new_track(1, name="Annotated Features", scale=0)
            feature_set = track_for_features.new_set()
        else:
            raise ValueError(f"You need to pass both a diagram and feature_set.")

        #start_feature_list = []
        #end_feature_list = []
        for locustag, feature in self:
            feat_sigil = "ARROW"
            feat_arrowshaft = 0.6
            feat_border_color = colors.black
            if key_color:
                feat_color = key_color(*key_color_args)
            else:
                feat_color = colors.white

            #start_feature_list.append(feature.start)
            #end_feature_list.append(feature.end)

            feature_set.add_feature(feature.feature, sigil=feat_sigil, color=feat_color, border=feat_border_color,  arrowshaft_height=feat_arrowshaft)

        diagram.draw(format=format, pagesize=pagesize, fragments=1, fragment_size=0.99, track_size=0.6, start=self.left_origin, end=self.right_origin, x=0.01, y=0.01)
        diagram.write(f"{output}.{output_format}", output_format)

    def draw_EGF(self, output):
        self.reset_origin_coordinate()
        features=[]
        for locustag, feature in self:
            #Add color management
            features.append(GraphicFeature(start=feature.start, end=feature.end, strand=feature.strand, color="#ffcccc", label=locustag))

        record = GraphicRecord(first_index=self.left_origin, sequence_length=len(self.seq[0]), features=features)
        record.plot(figure_width=8)
        plt.savefig(output, dpi=300)



################# Functions

def process_genome(genbank, tag=None):
    if tag:
        return (tag, Genome(genbank))
    else:
        return Genome(genbank)

def process_multiple_genome(genbanks, tags=None, process=4):
    with Pool(processes=process) as pool:
        if tags and len(tags) == len(genbanks):
            results = [pool.apply_async(process_genome, args=(gbk, t)) for gbk, t in zip(genbanks, tags)]
            output = [p.get() for p in results]
        else:
            results = [pool.apply_async(process_genome, args=(gbk, )) for gbk in genbanks]
            output = [p.get() for p in results]
    
    return output

def parallel_genomes(genbank_files, workers=8, print_out=False):
    with ProcessPoolExecutor(max_workers=workers) as executor:
        # Iterator to remember where we are in the list of genbank
        gbks_iter = iter(genbank_files)
        # Initial submit
        futures = {executor.submit(process_genome, gbk) for gbk in list(itertools.islice(gbks_iter, workers))}

        # All done when no other job is submited
        while futures:
            # Wait for done processes and replace futures list without completed processes
            done, futures = concurrent.futures.wait(futures, return_when=concurrent.futures.FIRST_COMPLETED)

            # Resubmit len(done) job
            for gbk in list(itertools.islice(gbks_iter, len(done))):
                futures.add(executor.submit(process_genome, gbk))

            # process done job
            for future in done:
                genome = future.result()
                strain = genome.strain
                if print_out:
                    print(f"{strain}")
                yield (strain, genome)

def sequential_genomes(genbank_files):
    for gbk in genbank_files:
        genome = Genome(gbk)
        strain = genome.strain
        yield (strain, genome)
