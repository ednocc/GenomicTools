import os
import glob
import subprocess
import pandas

from genomictools.gff.GFFParser import GFFParser
from genomictools.database.Database import Database
from genomictools.hashfile.HashFile import HashFile

__all__ = ["Roary"]

class Roary:
    def __init__(self, roary_dir, gff_dir, database_path):
        self.roary_dir = os.path.realpath(roary_dir)
        self.roary_content = os.listdir(self.roary_dir)
        self.list_gff = glob.glob(f"{gff_dir}*.gff")
        self.renamed = False
        self.acc_to_strain = dict(Database.convert_accession_to_strain(database_path))

    def correct_roary_table(self):
        accessions = [GFFParser(gff).gff_name for gff in self.list_gff]
        fds = [open(gff, "rb") for gff in self.list_gff]
        hasher = HashFile()
        hashes = [(accession, 
                acc_to_strain[accession],
                gff,
                hasher.hash_file(fd)) for accession, gff, fd in zip(accessions, self.list_gff, fds)]
        [fd.close() for fd in fds]
        del fds
        #print(hashes)

        # Use accession or strain name in the pandas table
        if self.renamed:
            i = 1
        else:
            i = 0

        try:
            table = self.roary_content.index("gene_presence_absence.csv")
        except ValueError as e:
            raise FileNotFoundError(f"File 'gene_presence_absence.csv' not in directory {self.roary_dir}")
            
        table = f"{self.roary_dir}/gene_presence_absence.csv"
        csv = pandas.read_csv(table)
        for h in hashes:
            print("\t".join(h))
            with open(h[2], "r") as g:
                rows = list(GFFParser(g).parse_gff())
                csv[h[i]] = csv[h[i]].apply(Roary._apply_convert, args=[h, rows])

        csv.to_csv(f"{self.roary_dir}/gene_presence_absence_modified.csv", index=False, sep=",")
    
    # Takes on input each element of the column individualy and return the modified element
    @staticmethod
    def _apply_convert(elem, h, rows):
        try:
            if str(elem) != "nan":
                multi = elem.split("\t")
                temp = []
                for m in multi:
                    if m.split("_")[0] == h[3]:
                        # gff index start at 1 but list index start at 0 
                        # so we need to remove 1 to the gff index in order 
                        # to keep the correct list index
                        index = int(m.split("_")[1]) - 1
                        locus = rows[index][-1]["ID"]
                        temp.append(locus)
                elem = ";".join(temp)
        except Exception as e:
            print(e)
            print(elem)
            exit()
    
        return elem

    def translate_accession(self, file_to_translate):
        path = os.path.realpath(file_to_translate)
        s = os.path.splitext(path)
        exe = os.path.basename(s[0])
        extension = s[1]
        #if extension in [".newick", ".nwk"]:
        #    output = f"{self.roary_dir}/{file_to_translate.replace('.newick', '_translate.newick')}"
        #elif extension == ".csv":
        #    output = f"{self.roary_dir}/{file_to_translate.replace('.csv', '_translate.csv')}"
        #else:
        output = f"{self.roary_dir}/{exe} + '_translate' + {extension}"

        with open(path, "r") as fd:
            with open(output, "w") as out:
                temp = fd.read()
                for accession, strain in self.acc_to_strain.items():
                    #f_ext = os.path.splitext(str)
                    #name = "_".join(self.gbk.strip(f_ext).split("_")[:3])
                    if "_genomic.gbff" in temp:
                        temp = temp.replace(accession + "_genomic.gbff", strain)
                    elif "_true_genomic.gbk" in temp:
                        temp = temp.replace(accession + "_true_genomic.gbk", strain)
                    elif ".gbk" in temp:
                        temp = temp.replace(accession + ".gbk", strain)
                    else:
                        temp = temp.replace(accession, strain)
                out.write(temp)
    
    def roary_plot(self):
        exe_dir = os.path.split(__file__)[0]
        tree = f"{self.roary_dir}/accessory_binary_genes.fa"
        table = f"{self.roary_dir}/gene_presence_absence.csv"
        call = subprocess.call(["python3", f"{exe_dir}/roary_plots.py", "--labels", tree, table])

