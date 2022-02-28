import pandas
import numpy as np
import sys
import os

def main():
    output = sys.argv[1]
    vcf = sys.argv[2:]
    df = []
    for v in vcf:
        with open(v, 'r') as vcf_in, open(f"temp_{v}", "w") as out_vcf:
            for row in vcf_in.read().splitlines():
                start_line = row[:2]
                if start_line != "##":
                    print(row, file=out_vcf)
        df.append([v, pandas.read_csv(f"temp_{v}", sep="\t")])
        os.remove(f"temp_{v}")

    df_index = []
    named_index = []
    for v, d in df:
        d.set_index("POS", inplace=True, drop=True)
        df_index.append(set(d.index))
        named_index.append([v, set(d.index)])

    with open(output, "w") as result:
        intersection = set.intersection(*df_index)
        print("Intersection : ", len(intersection))
        print("Intersection : ", len(intersection), file=result)
        print(sorted(list(intersection)), file=result)

        for lv, l in named_index:
            for rv, r in named_index:
                if lv == rv:
                    continue
                else:
                    print(f"{lv} vs {rv} : {len(l - r )}")
                    print(f"{lv} vs {rv} : {len(l - r )}", file=result)
                    print(sorted(list(l - r)), file=result)

    #run1 = pandas.read_csv(sys.argv[1], sep="\t")
    #run1bis = pandas.read_csv(sys.argv[2], sep="\t")
    #merge = pandas.read_csv(sys.argv[3], sep="\t")

    #run1.set_index("POS", inplace=True, drop=True)
    #run1bis.set_index("POS", inplace=True, drop=True)
    #merge.set_index("POS", inplace=True, drop=True)
    #
    ##run1 = run1[["TYPE", "REF", "ALT", "EVIDENCE", "LOCUS_TAG"]]
    ##run1 = run1[run1["TYPE"] == "snp"]

    ##run1bis = run1bis[["TYPE", "REF", "ALT", "EVIDENCE", "LOCUS_TAG"]]
    ##run1bis = run1bis[run1bis["TYPE"] == "snp"]

    ##merge = merge[["TYPE", "REF", "ALT", "EVIDENCE", "LOCUS_TAG"]]
    ##merge = merge[merge["TYPE"] == "snp"]

    ##print("Run1")
    ##print(run1)
    ##print("Run1bis")
    ##print(run1bis)
    ##print("Merge")
    ##print(merge)

    #index_run1 = set(run1.index)
    #index_run1bis = set(run1bis.index)
    #index_merge = set(merge.index)

    #intersection = index_run1 & index_run1bis & index_merge
    #print(len(intersection))
    #print(sorted(list(intersection)))

    #print(len(index_run1 - index_run1bis))
    #print(len(index_run1 - index_merge))
    #print(len(index_run1bis - index_run1))
    #print(len(index_run1bis - index_merge))
    #print(len(index_merge - index_run1))
    #print(len(index_merge - index_run1bis))

if __name__ == "__main__":
    main()
