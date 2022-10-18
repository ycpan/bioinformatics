from BCBio import GFF
#in_handle = open("tunisia.gff",'r')
#in_handle = open("GCF_000001405.40_GRCh38.p14_genomic.gff",'r')
in_handle = open("ANXA1.gff",'r')
for rec in GFF.parse(in_handle):
    for feature in rec.features:
        if feature.type == "region":
            chromo = feature.qualifiers['ID'][0]
            chrid = chromo.split(":")[0]
            chrLen = chromo.split("..")[1]
            print(chrid,chrLen)
