import re
import sys

# Input SAM files from bowtie or bowtie2 to call SNP and its reads number;
# Demand the reads that have SNP excluding INDEL;

inputfile = sys.argv[1]     #input first arg: .sam file 
hash = {}

with open(inputfile, 'r') as f:
    for line in f:
        line = line.strip()
        match = re.search(r'MD:Z:((\d+[ATCGatcg])+)', line)
        if not match:
            continue
        snp_des = match.group(1)
        t = line.split('\t')
        cigar = t[5]
        seq = t[9]
        len_seq = len(seq)
        if str(len_seq) not in cigar:  # exclude reads with INDEL
            continue
        poffset = 0
        for m in re.finditer(r'(\d+)([ATCGatcg])', snp_des):
            offset = int(m.group(1)) + poffset
            snpRef = m.group(2)
            genomicpos = int(t[3]) + offset
            snpReads = seq[offset]
            poffset = offset + 1
            if snpReads == "N":
                continue
            des = snpRef + snpReads
            if (t[2], genomicpos) not in hash:
                hash[(t[2], genomicpos)] = {des: 1}
            else:
                if des not in hash[(t[2], genomicpos)]:
                    hash[(t[2], genomicpos)][des] = 1
                else:
                    hash[(t[2], genomicpos)][des] += 1

# input samtools depth file
hashC = {}
with open(sys.argv[2], 'r') as f:  # input second arg: A.seqdepth.txt file
    for line in f:
        chr, pos, depth = line.strip().split('\t')
        hashC[(chr, int(pos))] = int(depth)

filename = inputfile.replace('.sam', '')
outfile = filename + '_all_SNP.bed'
with open(outfile, 'w') as out:
    for (chr, genomicpos) in sorted(hash):
        for des in hash[(chr, genomicpos)]:
            if (chr, genomicpos) not in hashC:
                hashC[(chr, genomicpos)] = 'NA'
                out.write(f"{chr}\t{genomicpos-1}\t{genomicpos}\t{des}\t{hash[(chr, genomicpos)][des]}\t{hashC[(chr, genomicpos)]}\t'NA'\n")    #output .bed format for bedtools intersect
            else:
                percent = 100 * (hash[(chr, genomicpos)][des] / hashC[(chr, genomicpos)])
                out.write(f"{chr}\t{genomicpos-1}\t{genomicpos}\t{des}\t{hash[(chr, genomicpos)][des]}\t{hashC[(chr, genomicpos)]}\t{percent:.4f}\n")
