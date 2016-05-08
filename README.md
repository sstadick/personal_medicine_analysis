
Sources:
- Cosmic curatate gene list: http://cancer.sanger.ac.uk/cosmic/curation
- VCF File: ftp://ftp2.completegenomics.com/vcf_files/Build37_2.0.0/vcfBeta-NA12880-200-37-ASM.vcf
- Liftover chain file: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/
- cosmic.vcf: http://cancer.sanger.ac.uk/cosmic/download -> /files/grch38/cosmic/v76/VCF/CosmicCodingMuts.vcf.gz and /files/grch38/cosmic/v76/VCF/CosmicNonCodingVariants.vcf.gz then append together

Plan:
- Using the pullCosm.py script found in this directory, I pulled the top mutations for each fo the 178 genes in cosmic curated list. 
- iterate through vcf file and look for a mutation that matches one of the 178 hotspot mutations
- results can be found in ./mutationsfound.vcf

Read in the 178 genes (2 were lost in liftover)


```python
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
```

    ['chr9', '133589333', '133763062', 'ABL1']
    ['chr2', '158592956', '158732374', 'ACVR1']
    ['chr12', '52345451', '52390862', 'ACVR1B']
    ['chr14', '105235686', '105262088', 'AKT1']
    ['chr2', '29415640', '30144432', 'ALK']
    ['chrX', '63404997', '63425624', 'AMER1']
    ['chr5', '112043195', '112181936', 'APC']
    ['chr1', '27022524', '27108595', 'ARID1A']
    ['chr12', '46123448', '46301823', 'ARID2']
    ['chr20', '30946147', '31027122', 'ASXL1']


Read in and check the vcf file at the same time


```python
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
```


```python
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
```


```python
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
```

    ['3', '37053568', '.', 'A', 'G', '.', '.', 'NS=1;AN=2;AC=1;CGA_XR=dbsnp.89|rs1799977;CGA_FI=4292|NM_000249.3|MLH1|CDS|MISSENSE&4292|NM_001167617.1|MLH1|CDS|MISSENSE&4292|NM_001167618.1|MLH1|UTR5|UNKNOWN-INC&4292|NM_001167619.1|MLH1|UTR5|UNKNOWN-INC', 'GT:PS:FT:GQ:HQ:EHQ:CGA_CEHQ:GL:CGA_CEGL:DP:AD:CGA_RDP', '1/0:.:PASS:249:366,249:366,249:53,50:-366,0,-249:-53,0,-50:32:18,14:14']
    ['3', '47125385', '.', 'G', 'A', '.', '.', 'NS=1;AN=2;AC=2;CGA_XR=dbsnp.108|rs4082155;CGA_FI=29072|NM_014159.6|SETD2|CDS|MISSENSE', 'GT:PS:FT:GQ:HQ:EHQ:CGA_CEHQ:GL:CGA_CEGL:DP:AD:CGA_RDP', '1/1:.:PASS:150:1218,150:1218,150:49,52:-1218,-150,0:-52,-49,0:56:56,56:0']
    ['3', '52584431', '.', 'T', 'TCTC', '.', '.', 'NS=1;AN=2;AC=1;CGA_XR=dbsnp.126|rs34372721&dbsnp.130|rs71084187;CGA_FI=55193|NM_018165.4|PBRM1|DONOR|UNKNOWN-INC&55193|NM_018313.4|PBRM1|DONOR|UNKNOWN-INC&55193|NM_181042.3|PBRM1|DONOR|UNKNOWN-INC', 'GT:PS:FT:GQ:HQ:EHQ:CGA_CEHQ:GL:CGA_CEGL:DP:AD:CGA_RDP', '1/0:.:PASS:292:292,396:290,395:39,41:-292,0,-396:-39,0,-41:35:20,15:15']
    ['4', '55972974', '.', 'T', 'A', '.', '.', 'NS=1;AN=2;AC=1;CGA_XR=dbsnp.92|rs1870377;CGA_FI=3791|NM_002253.2|KDR|CDS|MISSENSE', 'GT:PS:FT:GQ:HQ:EHQ:CGA_CEHQ:GL:CGA_CEGL:DP:AD:CGA_RDP', '1/0:.:PASS:316:430,316:427,313:50,53:-430,0,-316:-50,0,-53:42:25,17:17']
    ['5', '35871190', '.', 'G', 'A', '.', '.', 'NS=1;AN=2;AC=1;CGA_XR=dbsnp.88|rs1494555;CGA_FI=3575|NM_002185.2|IL7R|CDS|MISSENSE;CGA_PFAM=PFAM|PF00041|FN3', 'GT:PS:FT:GQ:HQ:EHQ:CGA_CEHQ:GL:CGA_CEGL:DP:AD:CGA_RDP', '1|0:35870814:PASS:226:226,226:221,221:39,49:-226,0,-226:-39,0,-49:40:16,24:24']



```python
print len(mutsinroi)
print mutsinroi[1]

```

    17
    ['3', '47125385', '.', 'G', 'A', '.', '.', 'NS=1;AN=2;AC=2;CGA_XR=dbsnp.108|rs4082155;CGA_FI=29072|NM_014159.6|SETD2|CDS|MISSENSE', 'GT:PS:FT:GQ:HQ:EHQ:CGA_CEHQ:GL:CGA_CEGL:DP:AD:CGA_RDP', '1/1:.:PASS:150:1218,150:1218,150:49,52:-1218,-150,0:-52,-49,0:56:56,56:0']



```python
with open("./mutationsfound.vcf", 'w') as out:
    for mut in mutsinroi:
        out.write("\t".join(mut) + "\n")
```


```python

```
