import csv
import sys
import argparse
import sqlite3

def compare_pos_from_vcf_and_database(database, vcf_pos):
    conn = sqlite3.connect(database)
    request = "SELECT DISTINCT position FROM snp"
    cur = conn.cursor()

    cur.execute(request)
    db_pos = [pos[0] for pos in cur.fetchall()]
    
    # Compter les doublons
    print("Doublons")
    doublons_db = dict([(n, db_pos.count(n)) for n in set(db_pos)])
    doublons_vcf = dict([(n, vcf_pos.count(n)) for n in set(vcf_pos)])
    print("db_pos : ", len(doublons_db))
    for k, v in doublons_db.items():
        if v > 1:
            print(f"{k} : {v}")
    print("vcf_pos : ", len(doublons_vcf))
    for k, v in doublons_vcf.items():
        if v > 1:
            print(f"{k} : {v}")

    print("Diff")

    # Recherche des différences
    uniq_vcf = set(vcf_pos) - set(db_pos)
    uniq_db = set(db_pos) - set(vcf_pos)
    comm = set(vcf_pos) & set(db_pos)
    all = set(vcf_pos) | set(db_pos)
    #print(uniq_vcf)
    print(len(uniq_vcf))
    #print(uniq_db)
    print(len(uniq_db))
    print(len(comm))
    print(len(all))


def check_snp_count(vcf, strains, database):
    motif_p = ".:.:.:.:.:.:.:."
    motif_v = "."
    i = 0
    snps = {}
    vcf_pos = []
    
    with open(strains, "r") as fds:
        samples = fds.readlines()

    with open(vcf, "r") as fdv:
            vcf_line = csv.reader(fdv, delimiter="\t", quotechar='"')
            for row in vcf_line:
                if "#" not in row[0]:
                    if len(str(row[3])) == 1 and len(str(row[4])) == 1 and "TYPE=snp" in row[7]:
                        vcf_pos.append(int(row[1]))
                        i += 1
                        for n, strain in enumerate(samples, 9):
                            if row[n] == motif_p or row[n] == motif_v:
                                continue
                            else:
                                if strain in snps.keys():
                                    snps[strain] += 1
                                else:
                                    snps[strain] = 1

    print(f"Total snp : {i}\n")
    for k, v in snps.items():
        print(f"{k} : {v}")
    
    if database:
        compare_pos_from_vcf_and_database(database, vcf_pos)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check SNP number in vcf combine or not (per sample or in total")
    parser.add_argument("-v", "--vcf", type=str, help="VCF file", required=True)
    parser.add_argument("-s", "--strains", type=str, help="Strains list file", required=True)
    parser.add_argument("-d", "--database", type=str, help="SNP database file", default="")
    #parser.add_argument("-n", type=str, help="Index position of the strain in the vcf file (start from 9 to {strain number} - 1 in combine vcf", default=9)
    args = parser.parse_args()

    check_snp_count(args.vcf, args.strains, args.database)
