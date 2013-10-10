import sys, gzip, datetime
# Author: Liping Hou 
# Email:  houliping@gmail.com
# Oct 1st, 2013

if len(sys.argv)!=6:
    print("Usage: python3 prob2plink.py minimac.info minimac.prob.gz out.fam out.dose bad_imputed.snps")
    sys.exit(1)

print()
print("Analysis started: {}".format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
print()
print("Reading minimac info file from [ {} ]".format(sys.argv[1])) 

# Set Rsq threshold to 0.5
Rsq = 0.5

SNP = []
AL1 = []
AL2 = []
genotyped = 0
imputed = 0
badimputedSNP = []
with open(sys.argv[1]) as info:
    for infol in info:
        res = infol.split()
        snp = res[0]
        al1 = res[1]
        al2 = res[2]
        SNP.append(snp)
        AL1.append(al1)
        AL2.append(al2)
        if(res[7] == '-' and float(res[6]) < Rsq):
            badimputedSNP.append(snp)
        if(res[7] == '-'):
            imputed += 1
        else:
            genotyped += 1
        

nsnp = len(SNP) - 1
nbadsnp = len(badimputedSNP)

print("Found {} SNPs from the info file, including {} genotyped SNPs and {} imputed SNPs".format(nsnp, genotyped - 1, imputed))
print("Writing {} bad imputed (Rsq < 0.5) SNPs to [ {} ]".format(nbadsnp, sys.argv[5]))

with open(sys.argv[5], 'w') as snp:
    for i in badimputedSNP:
        snp.write("{}\n".format(i))

FID = []
IID = []
PROB ={}
for i in range(1,nsnp+1):
    PROB[SNP[i]] = []

print("Reading minimac prob file from [ {} ]".format(sys.argv[2]))

with gzip.open(sys.argv[2]) as prob:
    for line in prob:
        tem = line.decode('utf-8').split()
        iid = tem[0].split('->')[0]
        fid = tem[0].split('->')[1]
        FID.append(iid)
        IID.append(fid)
        for n in range(1,nsnp+1):
            PROB[SNP[n]].extend(tem[n*2:(n*2 + 2)])
       
nind = len(FID)

print("{} individuals read from [ {} ]".format(nind, sys.argv[2]))
print("Converting minimac prob file to plink dosage file ...")
print("Writing pedigree information to [ {} ]".format(sys.argv[3]))
print("Gender and phenotype were set to missing (-9) in the fam file")
print("Please update the father IDs and Mother IDs if you have nonfounders")
print("Writing converted plink dosage file to [ {} ]".format(sys.argv[4]))

with open(sys.argv[3], 'w') as fam:
    for id in range(0, nind):
        fam.write("{}\t{}\t0\t0\t-9\t-9\n".format(FID[id],IID[id]))
with open(sys.argv[4],'w') as dose:
    for j in range(1, nsnp+1):
        dose.write("{}\t{}\t{}\t{}\n".format(SNP[j], AL1[j], AL2[j], "\t".join(str(i) for i in PROB[SNP[j]])))
print()
print("Analysis finished: {}".format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
print()
#print to standout
#for j in range(1,nsnp+1):
#    print(SNP[j], AL1[j], AL2[j], end=" ")
#    for geno in PROB[SNP[j]]:
#        print(geno, end=" ")
#    print('\n', end="")
