from genomictools.wrapper import GC

__all__ = ["Metadata"]

class Metadata:
    def __init__(self, genome):
        self.genome = genome
        self.source = self.genome.records[0].features[0]

    def extract_species(self):
        return self.source.qualifiers.get("organism")[0]
    
    
    def extract_strain(self):
        return self.source.qualifiers.get("strain")[0].replace(" ", "-")
        
    
    def extract_sequence_id(self):
        i = []
        for record in self.genome.records:
            i.append(self.genome.records[0].id)
        return i
    
    def extract_plasmid(self):
        if len(self.genome.records) > 1:
            return [p.id for p in self.genome.records[1:]]
        else:
            return []
    
    def count_gene(self):
        i = [0 for n in range(0, len(self.genome.records))]
        for n in range(0, len(self.genome.records)):
            for locustag, feature in self.genome:
                if feature['contig_id'] == n:
                    i[n] += 1
        return i
    
    def count_tRNA(self):
        i = [0 for n in range(0, len(self.genome.records))]
        for n in range(0, len(self.genome.records)):
            for locustag, feature in self.genome:
                if feature['contig_id'] == n and feature['type'] == 'tRNA':
                    i[n] += 1
        return i
    
    def count_rRNA(self):
        i = [0 for n in range(0, len(self.genome.records))]
        for n in range(0, len(self.genome.records)):
            for locustag, feature in self.genome:
                if feature['contig_id'] == n and featrue['type'] == 'rRNA':
                    i[n] += 1
        return i
    
    def count_coding_and_pseudo(self):
        i = [[0, 0] for n in range(0, len(self.genome.records))]
        for n in range(0, len(self.genome.records)):
            for locustag, feature in self.genome:
                if feature['contig_id'] == n and feature['type'] == 'CDS' and not feature['pseudo']:
                    i[n][0] += 1
                elif feature['contig_id'] == n and feature['pseudo']:
                    i[n][1] += 1
                else:
                    continue
        return i
    
    def count_hypothetical(self):
        i = [0 for n in range(0, len(self.genome.records))]
        for n in range(0, len(self.genome.records)):
            for locustag, feature in self.genome:
                if feature['contig_id'] == n and feature['product'] == "hypothetical protein":
                    i[n] += 1
        return i
    
    #@wrapper.timer
    def calc_GC_content(self):
        # 6 --> arbitrary cutoff (plasmid in genbank of complete genome)
        # Draft genome currently have more than 6 contigs
        if len(self.genome.records) < 7:
            return round(GC(self.genome.seq[0]), 2)
        else:
            pGC = 0
            i = 0
            for s in self.genome.seq:
                pGC += GC(s)
                i += 1
            pGC = pGC // i
            return round(pGC, 2)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return f"Metadata(genome={self.genome.name})"
