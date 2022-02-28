import pysam
import sys

if __name__ == '__main__':
    vcf = pysam.VariantFile(sys.argv[1])
    for rec in vcf.fetch():
        print(rec.pos)
        #print(rec.info["AO"])
        for sample, variantRecordSample in rec.samples.items():
            depth = variantRecordSample["DP"]
            if not depth or depth < 10:
                print("\t", sample, depth)

