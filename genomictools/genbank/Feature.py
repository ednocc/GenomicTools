import sys

from Bio.Seq import Seq
# Biopython 1.76
try:
    from Bio import Alphabet
except ImportError:
    Alphabet = None
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

__all__ = ["create_feature", "Feature"]

def create_feature(start, end, seq="ATGC", strand=None, type=None, qualifiers=None):
    location = FeatureLocation(start, end, strand)
    seqfeature = SeqFeature(location=location, type=type, qualifiers=qualifiers)
    if Alphabet:
        seq = Seq(seq, Alphabet.generic_dna)
    else:
        seq = Seq(seq)
    return Feature(seqfeature, seq)

class Feature:
    def __init__(self, seqfeature, seq, record_id, feat_id):
        self.seqfeature = seqfeature
        self.seq = seq
        self.record_id = record_id
        self.feat_id = feat_id
        if self.seqfeature.qualifiers.get('pseudo', False):
            self.pseudo = True
        else:
            self.pseudo = False

    ## Properties
    @property
    def type(self):
        return self.seqfeature.type

    @property
    def strand(self):
        return self.seqfeature.strand
    
    @property
    def start(self):
        return self.seqfeature.location.start.position
    
    @property
    def end(self):
        return self.seqfeature.location.end.position
    
    @property
    def qualifiers(self):
        return self.seqfeature.qualifiers
    
    @property
    def gene(self):
        return self.seqfeature.qualifiers.get('gene', ['NA'])[0]
    
    @property
    def locustag(self):
        return self.seqfeature.qualifiers.get('locus_tag', ['NA'])[0]

    @property
    def product(self):
        return self.seqfeature.qualifiers.get("product", ['NA'])[0]

    @property
    def protein(self):
        if self.pseudo:
            return None
        else:
            if self.seqfeature.qualifiers.get('translation', False):
                if Alphabet:
                    return Seq(self.seqfeature.qualifiers.get('translation')[0], Alphabet.generic_protein)
                else:
                    return Seq(self.seqfeature.qualifiers.get('translation')[0])
            else:
                # Cas des rRNA, tRNA, ncRNA, tmRNA...
                return None

    ## End

    def add_other(self, other):
        try:
            getattr(self, "other")
        except AttributeError:
            self.other = []

        self.other.append(other)
    
    def _getNewSeqFeature(self, start, end, strand, type, qualifiers):
        loc = FeatureLocation(start, end, strand=strand)
        newFeature = SeqFeature(loc, type=type, qualifiers=qualifiers)
        return newFeature

    def __add__(self, int_value):
        newStart = self.start + int_value
        newEnd = self.end + int_value
        newFeature = self._getNewSeqFeature(newStart, newEnd, self.strand, self.type, self.qualifiers)
        return Feature(newFeature, self.seq, self.record_id, self.feat_id)

    def __radd__(self, int_value):
        return self.__add__(int_value)

    def __sub__(self, int_value):
        newStart = self.start - int_value
        newEnd = self.end - int_value
        newFeature = self._getNewSeqFeature(newStart, newEnd, self.strand, self.type, self.qualifiers)
        return Feature(newFeature, self.seq, self.record_id, self.feat_id)

    def __rsub__(self, int_value):
        return self.__sub__(int_value)

    # For compatibility
    def lshift(self, shift):
        return self.__sub__(shift)

    def reverse_complement(self, length_seq):
        if self.strand == 1:
            newStrand = -1
        else:
            newStrand = 1

        newEnd = length_seq - self.start + 1
        newStart = length_seq - self.end + 1

        newFeature = self._getNewSeqFeature(newStart, newEnd, newStrand, self.type, self.qualifiers)
        return Feature(newFeature, self.seq.reverse_complement(), self.record_id, self.feat_id)

    def to_seqrecord(self, type="dna", id="", description=""):
        if not description:
            description = self.product

        if not id:
            id = self.locustag
            
        if type == "dna": # All feature have DNA sequence
            return SeqRecord(self.seq, id=id, description=description)
        elif type == "protein": # Not all feature have protein sequence (pseudogene, tRNA, ...)
            protein = self.protein
            if protein:
                return SeqRecord(protein, id=id, description=description)
            else:
                return None

    def to_fasta(self, handler=sys.stdout, **kwargs):
        to_close = False
        #call to_seqrecord with args and kwargs and write to handler
        try:
            getattr(handler, "write")
        except AttributeError:
            to_close = True
            handler = open(handler, "w")
        
        record = self.to_seqrecord(**kwargs)
        if record:
            SeqIO.write(record, handler, "fasta")
        else:
            print(f"{self.locustag} is a pseudogene")

        if to_close:
            handler.close()

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return f"Feature(locustag={self.locustag}, type={self.type}, start={self.start}, end={self.end}, product={self.product})"

    def __getitem__(self, key):
        if isinstance(key, str):
            return getattr(self, key)

    def __contains__(self, value):
        return True if self.start <= value <= self.end else False

    def __len__(self):
        return len(self.seq)

    def __eq__(self, other):
        return (self.start == other.start and 
                self.end == other.end and 
                self.seq == other.seq and 
                self.record_id == other.record_id)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __ge__(self, other):
        return self.start >= other.start

    def __gt__(self, other):
        return self.start > other.start

    def __le__(self, other):
        return self.start <= other.start

    def __lt__(self, other):
        return self.start < other.start

