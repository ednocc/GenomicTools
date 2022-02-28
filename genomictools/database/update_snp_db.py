#try:
#    import fire
#except ImportError:
#    import sys
#    print("conda activate fire-0.2.1")
#    sys.exit()

try:
    from GenomicPackage.src import snp
except ImportError:
    import sys
    sys.path.insert(0, "/projet/isp/galaxy_data/tools/")
    from GenomicPackage.src import snp


import argparse
import csv
import sqlite3


def update_snp_db(database, vcf, dryrun):
    conn = sqlite3.connect(database)
    cur_select = conn.cursor()
    cur_delete = conn.cursor()
    request_select = "SELECT * FROM snp GROUP BY position ORDER BY position"
    request_delete = "DELETE FROM snp WHERE position = ?"

    snps_pos = []
    with open(vcf, "r") as fd:
        for row in snp.parse_snp_vcf(fd):
            # Position
            snps_pos.append(int(row[1]))

    i = 0
    cur_select.execute(request_select)
    for row in cur_select.fetchall():
        if row[3] not in snps_pos:
            if dryrun:
                print(row)
                i += 1
            else:
                cur_delete.execute(request_delete, str(row[3]))
    
    if dryrun:
        print(f"Total number of entries to delete : {i}")
            


if __name__ == "__main__":
    #fire.Fire()
    parser = argparse.ArgumentParser(description="Delete entries from database which are not present in the vcf file")
    parser.add_argument("-d", "--database", type=str, help="Database to cure", required=True)
    parser.add_argument("-v", "--vcf", type=str, help="VCF file used to cure the database", required=True)
    parser.add_argument("-n", "--dryrun", action="store_true", help="Print each entry which will be delete without deleted it")
    args = parser.parse_args()

    update_snp_db(args.database, args.vcf, args.dryrun)
