import os
import re

#from pathlib import Path
#from collections import OrderedDict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC

from reportlab.lib import colors
#from reportlab.lib.units import cm

#from Bio.Graphics import GenomeDiagram
#from Bio.Graphics.GenomeDiagram import CrossLink

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

from multiprocessing import Pool

import itertools
import concurrent.futures
#from concurrent.futures import ProcessPoolExecutor, as_completed

FORWARD_STRAND = 1
REVERSE_STRAND = -1

__all__ = ["FORWARD_STRAND", "REVERSE_STRAND", "Genome", "process_multiple_genome", "parallel_genomes", "sequential_genomes"]


class GenomeIterator:
    def __init__(self, genome):
        self.genome = genome
        self.record_id = 0
        self.n_features_record = 0
        self.current_locus = 0
        self.cursor = 0
        self.max = len(genome) # Total number of locustag in all records

    def __iter__(self):
        return self

    #d = 0: {"L1": [1, 2], "L2": [3, 4, 5]}, 1: {"L3": [1, 2], "L4": [3, 4, 5]}}
    def __next__(self):
        if self.cursor >= self.max: # End of the file - StopIteration
            raise StopIteration
        self.n_features_record = len(self.genome.index[self.record_id])
        locustag = list(self.genome.index[self.record_id].keys())[self.current_locus]
        feature = self.genome.access_locustag(locustag)
        if feature:
            self.cursor += 1
            self.current_locus += 1
            if self.current_locus >= self.n_features_record: # End of this record - Next record
                self.record_id += 1
                self.current_locus = 0
            return (locustag, feature)


class Genome:
    # Prevents the direct use of this class
    #def __new__(cls):
    #    if BaseGenome.__name__ == cls.__name__:
    #        raise ValueError("You try to create an instance of the BaseGenome class which is not recommanded. Use the Genome class instead.")
    #    return super().__new__(cls)

    def __init__(self, genbank, seqrecords=None, left_origin=0, right_origin=0):
        self.gbk = os.path.basename(genbank)
        self.path = os.path.realpath(genbank)
        self.fasta = None
        self.extension = os.path.splitext(self.path)[1]
        self.name = "_".join(self.gbk.strip(self.extension).split("_")[:3])
        self.types = ["CDS", "tRNA", "rRNA", "ncRNA", "tmRNA", "misc_feature", "repeat_region", "gene"]
        #self.major_type = ["CDS", "tRNA", "rRNA", "ncRNA", "tmRNA"]
        #self.minor_type = ["misc_feature", "gene"]
        self.index = {}
        self.strain_sep = "-"
        self.left_origin = left_origin
        self.right_origin = right_origin
        if seqrecords:
            self.records = seqrecords
        else:
            self.records = []
        self._dict_index()

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        #return f"{self.__class__.__name__}(path={self.path})"
        return f"{self.__class__.__name__}(strain={self.strain}, species={self.species}, path={self.path}, gc={self.calc_GC_content()})"

    def __iter__(self):
        return GenomeIterator(self)

    def __len__(self):
        length = 0
        for record_id, index in self.index.items():
            length += len(index.keys())
        return length

    def __eq__(self, other):
        return self.index == other.index

    def __bool__(self):
        for locustag_dict in self.index.values():
            if locustag_dict:
                return True
        #return True if self.index else False
        return False

    def __add__(self, other):
        if not isinstance(other, Genome):
            raise ValueError(f"other must be instance of Genome, found {type(other)}")
        if len(self.records) > 1:
            raise ValueError(f"To concatenate two Genome, there must be only one sequence (number of sequence {len(self.sequence())}).")

        # Construct new seq
        seq = self.sequence() + other.sequence()

        # Extract seqfeature from self and other
        seqfeatures = self.features()
        _, _, other_selected_features = other.select_features(0, len(other.sequence()), 0)
        seqfeatures += [(feature + len(self.sequence())).seqfeature for feature in other_selected_features]

        # Construct new seqrecord based on new_index and new_seq
        seqrecord = self.get_region_seqrecord(seq, 0, seqfeatures)

        # Return new Genome instance
        return Genome(self.gbk, seqrecord)

    def parse_records(self):
        for record in SeqIO.parse(self.path, "genbank"):
            #self.records.append(record)
            yield record

    def _dict_index(self):
        append_records = True
        if self.records: # Not empty = Genome with seqrecords
            seqrecords = self.records
            append_records = False
        else: # Empty = New genome
            seqrecords = self.parse_records()

        for (record_id, record) in enumerate(seqrecords):
            if append_records:
                self.records.append(record)
            self.index[record_id] = {}
            # Skip source
            for (index, feature) in enumerate(record.features[1:], 1):
                try:
                    locustag = feature.qualifiers.get("locus_tag")[0]
                except TypeError as e:
                    continue
                if not locustag in self.index[record_id]:
                    self.index[record_id][locustag] = [index]
                else:
                    self.index[record_id][locustag].append(index)

    def translate_records_name(self, name):
        if isinstance(name, str):
            for record_id, record in enumerate(self.records):
                if record.id == name or record.name == name:
                    return record_id
            else:
                raise Exception(f"Record name or id was not be found in this genome. {record}")
        elif isinstance(name, int):
            if name >= len(self.records):
                raise ValueError(f"Record value is too high. Given: {name} Max: {len(self.records)}")
            return int(name)
        else:
            raise ValueError(f"Record need to be of type int or str. Given: {type(name)}")

    # Iter throught the underlying data structure
    def iterfeatures(self, record_id=0, feat_type=None):
        record_id = self.translate_records_name(record_id)
        # Skip source feature
        for feat_id, seqfeature in enumerate(self.features(record_id)[1:], 1):
            if feat_type and seqfeature.type != feat_type:
                continue
            yield Feature(seqfeature, seqfeature.extract(self.sequence(record_id)), record_id, feat_id, len(self.sequence(record_id)))

    def sort_features(self, features):
        return [feat for feat in sorted(features, key=lambda f: f.start)]

    def select_features(self, left, right, record_id=0, strict=True, feat_type=None):
        left, right = int(left), int(right)
        selected_features = []
        # TODO: check if feature in selected_features
        # TODO: features with a CompoundLocation are always selected
        for feature in self.iterfeatures(record_id, feat_type):
            if left < feature.start < right:
                selected_features.append(feature)
            elif left < feature.end < right:
                selected_features.append(feature)
            else:
                continue

            #if l in feature and r in feature: # For region smaller than a feature
            #    selected_features.append(feature)
            #if feature.start > feature.end:
            #    if feature.start <= left < feature.end: # overlapping feature on the left border
            #else:
            #    if feature.start <= left < feature.end: # overlapping feature on the left border
            #        if not strict:
            #            left = feature.start
            #        selected_features.append(feature)
            #    elif left <= feature.start and feature.end <= right: # feature intern to the region
            #        selected_features.append(feature)
            #    elif feature.start < right and right <= feature.end:
            #        if not strict: # overlapping feature on the right border
            #            right = feature.end
            #        selected_features.append(feature)
            #    else:
            #        continue

        selected_features = self.sort_features(set(selected_features))
        #if not strict:
        #    new_start = selected_features[0].start
        #    if new_start > right:
        #        new_start = 0
        #    left = new_start
        #    new_end = selected_features[-1].end
        #    if new_end < :
        #        pass
        #    right = new_end

        return left, right, selected_features

    def access_seqfeature(self, left, right, record_id=0, strict=True, feat_type=None):
        left, right = int(left), int(right)

        if left >= right:
            right, _, selected_features_right = self.select_features(right, len(self.sequence(record_id)), record_id, strict, feat_type)
            # Substract left coordinate as seqfeatures start at 0.
            selected_features_right = [feat.lshift(left) for feat in selected_features_right]
            _, left, selected_features_left = self.select_features(0, left, record_id, strict, feat_type)
            # Add length from left to the end of the genome as the seqfeatures start directly after the right segment.
            selected_features_left = [(feat + (len(self.sequence(record_id)) - left)) for feat in selected_features_left]
            selected_features = selected_features_right + selected_features_left
        else:
            left, right, selected_features = self.select_features(left, right, record_id, strict, feat_type)
            selected_features = [feat.lshift(left) for feat in selected_features]

        # Sort features by starting position and convert to SeqFeature
        selected_features = [feat.seqfeature for feat in self.sort_features(set(selected_features))]

        return selected_features

    def get_seqfeature_by_rank(self, record_id, feat_id_list):
        types_rank = dict([(v, k) for k, v in enumerate(self.types)])
        seqfeatures_rank = []
        for feat_id in feat_id_list:
            seqfeature = self.features(record_id)[feat_id]
            rank = types_rank[seqfeature.type]
            seqfeatures_rank.append((feat_id, rank))
        # Sort list by rank and return the first feat_id
        selected_feat_id = sorted(seqfeatures_rank, key=lambda x: x[1])[0][0]
        seqfeature = self.features(record_id)[selected_feat_id]
        return selected_feat_id, seqfeature

    # Index based
    def access_locustag(self, locustag):
        for record_id in self.index:
            if self.index[record_id].get(locustag, False):
                feat_id_list = self.index[record_id][locustag]
                feat_id, seqfeature = self.get_seqfeature_by_rank(record_id, feat_id_list)
                return Feature(seqfeature, seqfeature.extract(self.sequence(record_id)), record_id, feat_id, len(self.sequence(record_id)))

    def access_feature_id(self, record_id, feat_id):
        record_id = self.translate_records_name(record_id)
        seqfeature = self.features(record_id)[feat_id]
        return Feature(seqfeature, seqfeature.extract(self.sequence(record_id)), record_id, feat_id, len(self.sequence(record_id)))

    def access_position(self, position, record_id=0):
        record_id = self.translate_records_name(record_id)

        try:
            for locustag, feature in self:
                if feature.record_id == record_id:
                    if position in feature:
                        # If feature len feature is the same as record_id len so this is probably an error
                        if len(feature) == len(self.records[record_id]):
                            continue
                        return feature
                    else:
                        if position < feature.start:
                            return "Intergenic"
            else: # Intergenic after the last gene
                return "Intergenic"
        except Exception as e:
            raise ValueError(f"record_id params must be int or str. {type(record_id)}.\n({e})")

    def access_seqrecord(self, record_id):
        record_id = self.translate_records_name(record_id)
        return Genome(self.gbk, [self.records[record_id]])

    def get_region_seqrecord(self, seq, record_id, seqfeatures):
        seqrecord = [SeqRecord(seq, 
                               id=self.records[record_id].id, 
                               name=self.records[record_id].name, 
                               description=self.records[record_id].description, 
                               dbxrefs=self.records[record_id].dbxrefs, 
                               features=seqfeatures, 
                               annotations=self.records[record_id].annotations, 
                               letter_annotations=self.records[record_id].letter_annotations)]
        return seqrecord

    def get_region_sequence(self, left, right, record_id=0):
        left, right = int(left), int(right)

        # If right coordinate is outside the sequence then assume that this mean user want the beginning of the circular chromosome
        #if right >= len(self.sequence(record_id)):
        #    right %= len(self.sequence(record_id))

        if left >= right:
            seq_right = self.sequence(record_id)[int(right):len(self.sequence(record_id))]
            seq_left = self.sequence(record_id)[0:int(left)]
            seq = seq_right + seq_left
        else:
            seq = self.sequence(record_id)[int(left):int(right)]

        if not len(seq):
            raise ValueError(f"Extracted sequence is empty using coordinates {left}:{right} on record_id {record_id} ({self.records[record_id].id}, {len(self.sequence(record_id))}).")

        return seq

    def reverse_complement_record(self, record_id=0):
        record_id = self.translate_records_name(record_id)

        new_seq = self.get_region_sequence(0, len(self.sequence(record_id)), record_id=record_id)
        new_seq = new_seq.reverse_complement()
        seqfeatures = [self.source(record_id)]
        seqfeatures += [feat.reverse_complement().seqfeature for feat in self.iterfeatures(record_id)]

        return self.get_region_seqrecord(new_seq, record_id, seqfeatures)

    def convert_coordinate_on_circular_genome(self, position, record_id=0):
        record_id = self.translate_records_name(record_id)

        if position < 0 or position >= len(self.sequence(record_id)):
            return position % len(self.sequence(record_id))
        else:
            return position

    def access_region(self, left, right, record_id=0, strict=True, reverse=False, feat_type=None):
        left, right = int(left), int(right)
        record_id = self.translate_records_name(record_id)

        if reverse:
            new_record_id = 0
            new_genome = Genome(self.gbk, self.reverse_complement_record(record_id))
            # TODO: not working
            new_left = len(new_genome.sequence(new_record_id)) - right
            new_right = len(new_genome.sequence(new_record_id)) - left
            reverse = False
            return new_genome.access_region(new_left, new_right, new_record_id, strict, reverse, feat_type)

        # Sequence
        new_seq = self.get_region_sequence(int(left), int(right), record_id=record_id)

        # Features
        # Add source feature
        seqfeatures = [self.source(record_id)]
        seqfeatures += self.access_seqfeature(left, right, record_id, strict, feat_type)

        # Reverse
        #if reverse:
        #    new_seq = new_seq.reverse_complement()
        #    seqfeatures += [feat.reverse_complement(len(new_seq)).seqfeature for feat in selected_features]
        #else:
        #    seqfeatures += [feat.seqfeature for feat in selected_features]

        return self.get_region_seqrecord(new_seq, record_id, seqfeatures)

    def extract_region(self, left, right, record_id=0, strict=True, reverse=False, feat_type=None):
        seqrecord = self.access_region(left, right, record_id, strict, reverse, feat_type)
        return Genome(self.gbk, seqrecord, left, right)

    def __getitem__(self, key):
        if isinstance(key, str): # Single locustag
            if not self:
                raise ValueError(f"Genome index is empty. {self.index}")
            return self.access_locustag(key)
        elif isinstance(key, slice) and isinstance(key.start, str) and isinstance(key.stop, str): # From locustag to locustag (l1:l2)
            feat1 = self.access_locustag(key.start)
            if not feat1:
                raise KeyError(f"{key.start} not in genome.")

            feat2 = self.access_locustag(key.stop)
            if not feat2:
                raise KeyError(f"{key.stop} not in genome.")

            if feat1.record_id != feat2.record_id:
                raise KeyError(f"Passed keys are not on the same sequence. Start: {feat1.locustag} {feat1.record_id}. Stop : {feat2.locustag} {feat2.record_id}")
            else:
                record_id = feat1.record_id

            start, end = feat1.start, feat2.end
        elif isinstance(key, tuple):
            if isinstance(key[1], int): # (record_id, position)
                record_id, position = key
                return self.access_position(position, record_id)
            elif isinstance(key[1], slice): # From position to position (record_id, p1:p2)
                # No need to check for a particular type as translate_records_name is called in access_region
                record_id = key[0]
                start = key[1].start
                end = key[1].stop
                if not (isinstance(start, int) and isinstance(end, int)):
                    raise KeyError(key)
            else:
                raise KeyError(key)
        else:
            raise KeyError(key)

        #if start > end:
        #    start, end = end, stop
        seqrecord = self.access_region(start, end, record_id)

        return Genome(self.gbk, seqrecord, start, end)

    ######  Search  #######

    def search_gene(self, gene):
        for locustag, feature in self:
            if gene == feature.gene:
                return feature

    def search_product(self, product, regex=False):
        for locustag, feature in self:
            if regex:
                s = re.search(product, feature.product)
                if s:
                    return feature
            else:
                if product in feature.product:
                    return feature

    def extract_locustag_before_after(self, start, end, record_id=0, pseudo=True):
        before = []
        after = []
        record_id = self.translate_records_name(record_id)

        for locustag, feature in self:
            if feature.record_id == record_id:
                if feature.start < start:
                    before.append((locustag, feature))
                if feature.start > end:
                    after.append((locustag, feature))

        # Filtering
        if not pseudo:
            before = [(l, f) for l, f in before if not f.pseudo]
            after = [(l, f) for l, f in after if not f.pseudo]

        ret = dict([before[-1], after[0]])
        return ret

    def change_origin(self, record_id=0, gene_name=None, locustag=None, position=None, reverse=False, debug=False):
        # Find the locustag or gene_name
        if gene_name:
            starting_feature = self.search_gene(gene_name)
            if starting_feature:
                if reverse:
                    start = starting_feature.end
                else:
                    start = starting_feature.start
            else:
                raise ValueError(f"Gene {gene_name} not found. Can not change origin.")
        if locustag:
            starting_feature = self[locustag]
            if reverse:
                start = starting_feature.end
            else:
                start = starting_feature.start
        if position:
            start = position

        #seg_right = self.extract_region(start, len(self.sequence(record_id)), record_id=record_id, reverse=reverse)
        #seg_left = self.extract_region(0, start, record_id=record_id, reverse=reverse)
        new_genome = self.extract_region(start, start, record_id=record_id, reverse=reverse)

        #Debug
        #if debug:
        #    seg_left.format("Debug_left.gbk", "genbank")
        #    print("Left segment features:")
        #    for feature in seg_left.iterfeatures():
        #        print(feature)
        #    seg_right.format("Debug_right.gbk", "genbank")
        #    print("Right segment features:")
        #    for feature in seg_right.iterfeatures():
        #        print(feature)

        #new_genome = seg_right + seg_left
        return new_genome

    #####  Drawing  #####

    def reset_origin_coordinate(self, record_id=0):
        record_id = self.translate_records_name(record_id)
        self.records[record_id].features = [self.source(record_id)] + [(feature + self.left_origin).seqfeature for feature in self.iterfeatures(record_id)]
        #for feature in self.iterfeatures():
        #    feature = feature + self.left_origin

    def apply_new_coordinate(self, record_id=0):
        record_id = self.translate_records_name(record_id)
        self.records[record_id].features = [self.source(record_id)] + [(feature - self.left_origin).seqfeature for feature in self.iterfeatures(record_id)]
        #for feature in self.iterfeatures():
        #    feature = feature + self.left_origin

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

        record = GraphicRecord(first_index=self.left_origin, sequence_length=len(self.sequence()), features=features)
        record.plot(figure_width=8)
        plt.savefig(output, dpi=300)

    # Generic function
    def format(self, output, fmt, record_id=None):
        if record_id:
            record_id = self.translate_records_name(record_id)
            records = self.records[record_id]
        else:
            records = self.records
        count = SeqIO.write(records, output, fmt)
        return output

    # Predefined format
    def to_fasta(self, output, record_id=None):
        return self.format(output, "fasta", record_id=record_id)

    def to_genbank(self, output, record_id=None):
        return self.format(output, "genbank", record_id=record_id)

    def check_fasta(self, output=""):
        if not output:
            self.fasta = self.path.replace(self.extension, '.fna')
        else:
            self.fasta = output

        if not os.path.exists(self.fasta):
            self.to_fasta(self.fasta)

        return self.fasta

    #####  Metadata  #####

    def source(self, record_id=0):
        record_id = self.translate_records_name(record_id)
        return self.records[record_id].features[0]

    def features(self, record_id=0):
        record_id = self.translate_records_name(record_id)
        return self.records[record_id].features

    def sequence(self, record_id=0):
        record_id = self.translate_records_name(record_id)
        return self.records[record_id].seq

    def dbxref(self, record_id=0):
        record_id = self.translate_records_name(record_id)
        return self.records[record_id].features[0].qualifiers.get("db_xref", ["NA"])[0]

    @property
    def get_sequences_id(self):
        return [record.id for record in self.records]

    @property
    def species(self):
        species = self.source().qualifiers.get("organism")[0]
        if "subsp." in species:
            return " ".join(species.split()[:4])
        else:
            return " ".join(species.split()[:2])

    @property
    def strain(self):
        name_fields = ["strain", "isolate"]
        for field in name_fields:
            strain_name = self.source().qualifiers.get(field, None)
            if strain_name:
                return strain_name[0].replace(" ", self.strain_sep).replace("/", self.strain_sep)

        # If no strain name in source then return first record id
        return self.get_sequences_id[0]
        #raise ValueError("No field name found.")
        #return self.source.qualifiers.get("strain")[0].replace(" ", self.strain_sep).replace("/", self.strain_sep)

    @property
    def short_strain(self):
        #return self.source.qualifiers.get("strain")[0].replace(" ", "").replace("/", "").replace("-", "").replace("_", "")
        return self.strain.replace("-", "").replace("_", "")

    @property
    def species_code(self):
        return self.species.split()[-1][:4]

    def calc_GC_content(self):
        # 6 --> arbitrary cutoff (plasmid in genbank of complete genome)
        # Draft genome currently have more than 6 records
        if len(self.records) < 7:
            return round(GC(self.sequence()), 2)
        else:
            pGC = 0
            i = 0
            for r in self.records:
                pGC += GC(r.seq)
                i += 1
            pGC = pGC // i
            return round(pGC, 2)

    #def extract_plasmid(self):
    #    if len(self.records) > 1:
    #        return [p.id for p in self.records[1:]]
    #    elif len(self.records) > 20:
    #        raise ValueError("Your genome seems to be an assembly. This method can be use only on complete genome.")
    #    else:
    #        return []
    #
    #def count_gene(self):
    #    i = [0 for n in range(0, len(self.records))]
    #    for n in range(0, len(self.records)):
    #        for locustag, feature in self:
    #            if feature.record_id == n:
    #                i[n] += 1
    #    return i
    #
    #def count_tRNA(self):
    #    i = [0 for n in range(0, len(self.records))]
    #    for n in range(0, len(self.records)):
    #        for locustag, feature in self:
    #            if feature.record_id == n and feature.type == 'tRNA':
    #                i[n] += 1
    #    return i
    #
    #def count_rRNA(self):
    #    i = [0 for n in range(0, len(self.records))]
    #    for n in range(0, len(self.records)):
    #        for locustag, feature in self:
    #            if feature.record_id == n and feature.type == 'rRNA':
    #                i[n] += 1
    #    return i
    #
    #def count_coding_and_pseudo(self):
    #    i = [[0, 0] for n in range(0, len(self.records))]
    #    for n in range(0, len(self.records)):
    #        for locustag, feature in self:
    #            if feature.record_id == n and feature.type == 'CDS' and not feature.pseudo:
    #                i[n][0] += 1
    #            elif feature.record_id == n and feature.pseudo:
    #                i[n][1] += 1
    #            else:
    #                continue
    #    return i
    #
    #def count_hypothetical(self):
    #    i = [0 for n in range(0, len(self.records))]
    #    for n in range(0, len(self.records)):
    #        for locustag, feature in self:
    #            if feature.record_id == n and feature.product == "hypothetical protein":
    #                i[n] += 1
    #    return i
    

#class Genome(BaseGenome):
#    #@timer
#    def __init__(self, genbank):
#        super().__init__(genbank)
#        self._dict_index()
#
#    def __repr__(self):
#        return f"{self.__class__.__name__}(strain={self.strain}, species={self.species}, path={self.path}, gc={self.calc_GC_content()})"
#
#    ###### Search #######
#
#    def extract_region(self, left, right, record_id=0, strict=True, reverse=False):
#        seqrecord = self.access_region(left, right, record_id, strict, reverse)
#        return PartialGenome(self.gbk, seqrecord, left, right)
#
#    def dnaA_at_origin(self, output=""):
#        if output:
#            new_fna = Path(self.check_fasta(output))
#        else:
#            new_fna = Path(self.check_fasta() + ".tmp")
#
#        dnaA_locustag, dnaA_feature = [(locustag, feature) for locustag, feature in self.search_product("DnaA")]
#        if not type(dnaA_feature) == list:
#            dnaA_start = dnaA_feature.start
#            if dnaA_start not in range(0, 11):
#                self.seq[0] = self.seq[0][dnaA_start:] + self.seq[0][:dnaA_start]
#                self.to_fasta(new_fna)
#                return new_fna
#            else:
#                return self.fasta
#        else:
#            raise TypeError(f"'dnaA_feature' variable don't have the right type : {type(dnaA_feature)}")
#        
#
#    def change_origin(self, locustag=None, gene_name=None):
#        # Find the locustag or gene_name
#        if gene_name:
#            starting_feature = self.search_gene(gene_name)
#            start = starting_feature.start
#        if locustag:
#            starting_feature = self[locustag]
#            start = starting_feature.start
#
#        seg_right = self.extract_region(start, len(self.seq[0]))
#        seg_left = self.extract_region(0, start)
#        #Debug
#        #seg_left.format("Test_left.gbk", "genbank")
#        #seg_right.format("Test_right.gbk", "genbank")
#        new_genome = seg_right + seg_left
#        return new_genome
#
#
#    #def check_cached_data(self):
#    #    md5hash = self.md5hasher.hash_file(self.gbk)
#    #    path_cache_dir = Path(__file__).resolve().parent.joinpath(".cache")
#    #    if not path_cache_dir.is_dir():
#    #        path_cache_dir.mkdir()
#    #    cache_list = [f.name for f in path_cache_dir.glob("*")]
#    #    path_cache_file = path_cache_dir.joinpath(md5hash)
#    #    if md5hash in cache_list:
#    #        self.cached = True
#    #        with open(path_cache_file, "rb") as fd:
#    #            self.index = pickle.load(fd)
#    #    else:
#    #        # Save index object in binary file named with the md5 hash of the genbank file used to created it
#    #        self._dict_index()
#    #        with open(path_cache_file, "wb") as fd:
#    #            pickle.dump(self.index, fd)


#class PartialGenome(BaseGenome):
#    def __init__(self, genbank, seqrecord, left_origin=0, right_origin=0):
#        super().__init__(genbank)
#        self.left_origin = left_origin
#        self.right_origin = right_origin
#        self.records = seqrecord
#        self._dict_index()
#
#    def __add__(self, other):
#        if not isinstance(other, PartialGenome):
#            raise ValueError(f"other must be instance of PartialGenome, found {type(other)}")
#        if len(self.records) > 1:
#            raise ValueError(f"To concatenate two PartialGenome, there must be only one sequence (number of sequence {len(self.records[0].seq)}).")
#
#        # Construct new seq
#        seq = self.records[0].seq + other.records[0].seq
#
#        # Extract seqfeature from self and other
#        seqfeatures = self.records[0].features 
#        _, _, other_selected_features = other.selected_features(0, len(self.records[0].seq), 0)
#        seqfeatures += [(feature + len(self.records[0].seq)).seqfeature for feature in other_selected_features]
#
#        # Construct new seqrecord based on new_index and new_seq
#        seqrecord = self.get_region_seqrecord(seq, 0, seqfeatures)
#
#        # Return new PartialGenome instance
#        return PartialGenome(self.gbk, seqrecord)
#
#
#    def reset_origin_coordinate(self):
#        for locustag, feature in self:
#            feature = feature + self.left_origin
#
#    def draw(self, output, output_format, diagram=None, feature_set=None, diagram_name="MyDiagram", format="linear", pagesize=(1920, 80), key_color=None, key_color_args=None):
#        self.reset_origin_coordinate()
#
#        if not diagram and not feature_set:
#            diagram = GenomeDiagram.Diagram(diagram_name)
#            track_for_features = diagram.new_track(1, name="Annotated Features", scale=0)
#            feature_set = track_for_features.new_set()
#        else:
#            raise ValueError(f"You need to pass both a diagram and feature_set.")
#
#        #start_feature_list = []
#        #end_feature_list = []
#        for locustag, feature in self:
#            feat_sigil = "ARROW"
#            feat_arrowshaft = 0.6
#            feat_border_color = colors.black
#            if key_color:
#                feat_color = key_color(*key_color_args)
#            else:
#                feat_color = colors.white
#
#            #start_feature_list.append(feature.start)
#            #end_feature_list.append(feature.end)
#
#            feature_set.add_feature(feature.feature, sigil=feat_sigil, color=feat_color, border=feat_border_color,  arrowshaft_height=feat_arrowshaft)
#
#        diagram.draw(format=format, pagesize=pagesize, fragments=1, fragment_size=0.99, track_size=0.6, start=self.left_origin, end=self.right_origin, x=0.01, y=0.01)
#        diagram.write(f"{output}.{output_format}", output_format)
#
#    def draw_EGF(self, output):
#        self.reset_origin_coordinate()
#        features=[]
#        for locustag, feature in self:
#            #Add color management
#            features.append(GraphicFeature(start=feature.start, end=feature.end, strand=feature.strand, color="#ffcccc", label=locustag))
#
#        record = GraphicRecord(first_index=self.left_origin, sequence_length=len(self.seq[0]), features=features)
#        record.plot(figure_width=8)
#        plt.savefig(output, dpi=300)



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
