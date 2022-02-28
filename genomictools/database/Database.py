import glob
import csv
import os
import sqlite3
from pathlib import Path

from genomictools.wrapper import timer, pprint
from genomictools.genbank.Genome import Genome

__all__ = ["Database"]

class TableCreationError(Exception):
    """Exception use when tables in database are not create in the proper order"""

class Database:
    def __init__(self, path):
        self.path = Path()
        self.conn = sqlite3.connect(path)
        self.cursor = self.conn.cursor()
        self.path_to_strains_db = []

    # Close connection on program ending or when all reference to the object are delete
    def __del__(self):
        self.close()

    def close(self):
        self.cursor.close()
        self.conn.close()

    ########### Utility ###########

    def check_database(self):
        pass

    #def exists_table(self, table):
    #    check = f"SELECT name FROM slite_master WHERE type='table' AND name={table}"
    #    result = self.cursor.execute(check)
    #    return True if result.fetchall() else False

    def list_table(self):
        table = f"SELECT name FROM slite_master WHERE type='table'"
        result = self.cursor.execute(table)
        return result

    def list_path_to_strains_db(self, common_path):
        self.path_to_strains_db = glob.glob(f"{common_path}/FreeBayes/*.annoted.db")
   
    def convert_accession_to_strain(self):
        request = "SELECT accession, strain FROM strains"
        result = self.cur.execute(request)
        return [(accession, strain) for accession, strain in result.fetchall()]
    
    ########### End ###########

    ########### Table Creation ###########

    def create_strains_table(self):
        table_strains = """CREATE TABLE IF NOT EXISTS strains(
            id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
            strain TEXT UNIQUE NOT NULL,
            accession TEXT UNIQUE,
            species TEXT NOT NULL""" 

        self.cursor.execute(table_strains)
        self.conn.commit()

    def create_genbank_table(self):
        if not exists_table("strains"):
            raise TableCreationError("The strains table in the database must be create before genbank table.")

        table_genbank = """CREATE TABLE IF NOT EXISTS genbank(
            id INTEGER PROMARY KEY AUTOINCREMENT NOT NULL,
            key_strain TEXT NOT NULL,
            genome_size INTEGER NOT NULL,
            GCcontent REAL NOT NULL,
            uniq_gene INTEGER NOT NULL,
            coding_gene INTEGER NOT NULL,
            tRNA INTERGER NOT NULL,
            rRNA INTEGER NOT NULL,
            pre_locustag TEXT,
            path_gbk TEXT,
            FOREIGN KEY(key_strain) REFERENCES strains(strain) ON DELETE RESTRICT ON UPDATE CASCADE);"""

        self.cursor.execute(table_genbank)
        self.conn.commit()

    def create_refgenes_table(self):
        if not exists_table("strains"):
            raise TableCreationError("The strains table in the database must be create before refgenes table.")

        table_refgenes = """CREATE TABLE IF NOT EXISTS refgenes(
            id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
            key_strain TEXT NOT NULL,
            locustag INTEGER NOT NULL,
            type TEXT NOT NULL,
            start INTEGER NOT NULL,
            stop INTEGER NOT NULL,
            size INTEGER NOT NULL,
            pGC REAL,
            name TEXT,
            pseudo INTEGER,
            FOREIGN KEY(key_strain) REFERENCES strains(accession) ON DELETE RESTRICT ON UPDATE CASCADE);""" 
    
        self.cursor.execute(table_refgenes)
        self.conn.commit()


    def create_snp_table(self):
        if not exists_table("strains"):
            raise TableCreationError("The strains table in the database must be create before snp table.")

        if not exists_table("refgenes"):
            raise TableCreationError("The refgenes table in the database must be create before snp table.")

        table_snp = """CREATE TABLE IF NOT EXISTS snp(
            id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
            key_strain TEXT NOT NULL,
            chr TEXT NOT NULL,
            position INTEGER NOT NULL,
            ref TEXT NOT NULL,
            alt TEXT NOT NULL,
            qual REAL NOT NULL,
            key_locustag TEXT,
            annotation TEXT,
            effect TEXT,
            FOREIGN KEY(key_strain) REFERENCES strains(accession) ON DELETE RESTRICT ON UPDATE CASCADE,
            FOREIGN KEY(key_locustag) REFERENCES refgenes(locustag) ON DELETE RESTRICT ON UPDATE CASCADE);"""

        self.cursor.execute(table_snp)
        self.conn.commit()
    
    ########### End ###########

    ########### Data Insertion ###########

    def insert_into_snp_table(self, path_vcf_db, output):
        conn_vcf_db = sqlite3.connect(path_vcf_db)
        cursor_vcf_db = conn_vcf_db.cursor()
    
        strain_name = os.path.basename(path_vcf_db).replace(".final.annot.db", "")
    
        conn_snp_db = sqlite3.connect(output)
        cursor_snp_db = conn_snp_db.cursor()
    
        request_insert_into_snp_table = "INSERT INTO snp VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
        
        cursor_vcf_db.execute("SELECT CHROM, POS, REF, ALT, QUAL, ANN FROM alleles")
        for i in cursor_vcf_db.fetchall():
            if len(i[2]) == 1 and len(i[3]) == 1:
                if None in i:
                    data = (None, strain_name, *i[:len(i) - 1], None, None, None)
                    cursor_snp_db.execute(request_insert_into_snp_table, data)
                else:
                    annot = i[-1].split("|")
                    eff = ", ".join(annot[9:14])
                    data = (None, strain_name, *i[:len(i) - 1], annot[4], annot[1], eff)
                    cursor_snp_db.execute(request_insert_into_snp_table, data)
    
        conn_snp_db.commit()
        conn_snp_db.close()
        conn_vcf_db.close()
    
    
    def insert_into_refgenes_table(GBIndex, output):
        conn_snp_db = sqlite3.connect(output)
        cursor_snp_db = conn_snp_db.cursor()
        
        refname = os.path.basename(GBIndex.path).replace("_genomic.gbff", "")
    
        request_insert_into_refgenes_table = "INSERT INTO refgenes VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
    
        for k, v in GBIndex.index.items():
            for k2, v2 in v.items():
                type = v2.data.type
                start = v2.data.location.start
                stop = v2.data.location.end
                size = stop - start
                pGC = round(GC(v2.data.extract(GBIndex.records[0].seq)), 2)
                if v2.data.qualifiers.get('gene', []):
                    name = v2.data.qualifiers.get('gene')[0]
                else:
                    name = None
                if genbank.q.pseudo in v2.data.qualifiers:
                    pseudo = 1
                else:
                    pseudo = 0
                data = (None, refname, k, type, start, stop, size, pGC, name, pseudo)
                cursor_snp_db.execute(request_insert_into_refgenes_table, data)
    
        conn_snp_db.commit()
        conn_snp_db.close()
    
    
    def insert_into_strains_table(path_to_table, output):
        conn_snp_db = sqlite3.connect(output)
        cursor_snp_db = conn_snp_db.cursor()
    
        request_insert_into_strains_table = "INSERT INTO strains VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
    
        with open(path_to_table, "r") as table:
            reader = csv.reader(table, delimiter='\t', quotechar='"')
            for row in reader:
                data = (None, *row, None, None)
                cursor_snp_db.execute(request_insert_into_strains_table, data)
    
        conn_snp_db.commit()
        conn_snp_db.close()

    ########### End ###########
    
    
    #def convert_merge_vcf_to_sqlite3(merge_vcf):
    #    path = os.path.join(os.path.dirname(os.path.realpath(merge_vcf)), 'vcf.db')
    #    conn = sqlite3.connect("vcf.db")
    #
    #    if not os.path.exists(path):
    #        create_table = """create table vcf"""
    #        conn.execute(create_table)
    #
    #    with open(merge_vcf, "r") as fd:
    #        vcf_line = csv.reader(fd, delimiter='\t', quotechar='"')
    #        for row in vcf_line:
    #            if "#" not in row[0]:
    #                row = [str(x) for x in row]
    #                add_info = "insert into vcf values(" +  ", ".join(row) + ")"
    #                conn.execute(add_info)
    #
    #    conn.commit()
    #    conn.close() 
