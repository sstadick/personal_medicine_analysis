
# coding: utf-8

# Sources:
# - Cosmic curatate gene list: http://cancer.sanger.ac.uk/cosmic/curation
# - VCF File: ftp://ftp2.completegenomics.com/vcf_files/Build37_2.0.0/vcfBeta-NA12880-200-37-ASM.vcf
# - Liftover chain file: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/
# - cosmic.vcf: http://cancer.sanger.ac.uk/cosmic/download -> /files/grch38/cosmic/v76/VCF/CosmicCodingMuts.vcf.gz and /files/grch38/cosmic/v76/VCF/CosmicNonCodingVariants.vcf.gz then append together
# 
# Plan:
# - Using the pullCosm.py script found in this directory, I pulled the top mutations for each fo the 178 genes in cosmic curated list. 
# - iterate through vcf file and look for a mutation that matches one of the 178 hotspot mutations
# - results can be found in ./mutationsfound.vcf

# Read in the 178 genes (2 were lost in liftover)

# In[52]:

import sys
roiranges = []
with open("./ROI19.bed", 'r') as roi:
    for line in roi:
        if "#" not in line:
            elements = line.rstrip().split("\t")
            # chr, start, end
            elrange = [elements[0], elements[1], elements[2], elements[3]]
            roiranges.append(elrange)
            
for i in range(10):
    print roiranges[i]


# Read in and check the vcf file at the same time

# In[53]:

snps = []
indels = []
with open("./snp.bed", 'r') as snp:
    for line in snp:
        entry = line.rstrip().split("\t")
        snps.append(entry)

with open("./indel.bed", 'r') as indel:
    for line in indel:
        entry = line.rstrip().split("\t")
        indels.append(entry)


# In[54]:

def check_call_hot(call):
    for snp in snps:
        if "chr" + call[0] == snp[0]:
            if call[1] == snp[1]:
                return True
    for indel in indels:
        if "chr" + call[0] == indel[0]:
            if call[1] == indel[1]:
                return True
    return False


# In[55]:

mutsinroi = []
count = 0
with open("./vcfBeta-NA12880-200-37-ASM.vcf", 'r') as vcf:
    for line in vcf:
        if "#" not in line:
            count += 1
            call = line.rstrip().split("\t")
            if check_call_hot(call):
                mutsinroi.append(call)
                

for i in range(5):
    print mutsinroi[i]


# In[56]:

print len(mutsinroi)
print mutsinroi[1]


# In[59]:

with open("./mutationsfound.vcf", 'w') as out:
    for mut in mutsinroi:
        out.write("\t".join(mut) + "\n")


# In[ ]:



