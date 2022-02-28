from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

def add_annotation(seqrecord, start, end, strand=None, type="gene", qualifiers=None):
    feature = SeqFeature(location=FeatureLocation(start, end, strand=strand), type=type, qualifiers=qualifiers)
    seqrecord.features.append(feature)
    seqrecord.features = sorted(seqrecord.features, key=lambda feature: feature.location.start)
    return seqrecord
