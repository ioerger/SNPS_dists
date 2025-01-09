import sys

def read_genome(filename):
  s = ""
  n = 0
  for line in open(filename):
    if n==0: n = 1 # skip first
    else: s += line[:-1]
  return s

# 0-based
# start,end,Rv,gene,dir,descr

def read_genes(fname,descriptions=False):
  genes = []
  for line in open(fname):
    w = line.rstrip().split('\t')
    data = [int(w[1])-1,int(w[2])-1,w[8],w[7],w[3]]
    if descriptions==True: data.append(w[0])
    genes.append(data)
  return genes

def hash_genes(genes,genome):
  hash = {}
  for gene in genes:
    a,b = gene[0],gene[1]
    for i in range(a,b+1):
      hash[i] = gene
  return hash

def lookup_gene(i,genes_hash):
  return genes_hash.get(i,None)

codon    = { 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L',
             'CTA': 'L', 'CTG': 'L', 'CTN': 'L', 'TGG': 'W',
             'TAA': '*', 'TAG': '*', 'TGA': '*', 'ATG': 'M',
             'TTT': 'F', 'TTC': 'F', 'TAT': 'Y', 'TAC': 'Y',
             'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 
             'TCN': 'S', 'AGT': 'S', 'AGC': 'S', 'CCT': 'P', 
             'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CCN': 'P',
             'TGT': 'C', 'TGC': 'C', 'CAT': 'H', 'CAC': 'H',
             'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N',
             'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 
             'CGN': 'R', 'AGA': 'R', 'AGG': 'R', 'ATT': 'I', 
             'ATC': 'I', 'ATA': 'I', 'AAA': 'K', 'AAG': 'K',
             'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
             'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
             'ACN': 'T', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V',
             'GTG': 'V', 'GTN': 'V', 'GCT': 'A', 'GCC': 'A',
             'GCA': 'A', 'GCG': 'A', 'GCN': 'A', 'GGT': 'G', 
             'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'GGN': 'G',
             'TAN': 'X', 'TTN': 'X', 'TGN': 'X', 'CAN': 'X', 
             'ATN': 'X', 'AAN': 'X', 'GAN': 'X', 'AGN': 'X',
             'ANA': 'X', 'ANT': 'X', 'ANG': 'X', 'ANC': 'X', 
             'TNA': 'X', 'TNT': 'X', 'TNG': 'X', 'TNC': 'X', 
             'GNA': 'X', 'GNT': 'X', 'GNG': 'X', 'GNC': 'X', 
             'CNA': 'X', 'CNT': 'X', 'CNG': 'X', 'CNC': 'X', 
             'NAA': 'X', 'NAT': 'X', 'NAG': 'X', 'NAC': 'X', 
             'NTA': 'X', 'NTT': 'X', 'NTG': 'X', 'NTC': 'X', 
             'NGA': 'X', 'NGT': 'X', 'NGG': 'X', 'NGC': 'X', 
             'NCA': 'X', 'NCT': 'X', 'NCG': 'X', 'NCC': 'X', 
             'NNA': 'X', 'NNT': 'X', 'NNG': 'X', 'NNC': 'X', 
             'ANN': 'X', 'TNN': 'X', 'GNN': 'X', 'NNC': 'X', 
             'NAN': 'X', 'NTN': 'X', 'NGN': 'X', 'NCN': 'X', 
             'NNN': 'X'}

# genomepos and gene start/end (from annotation hash) are 0-based
# returns (Rv, genename, nucpos, tot, frame, refcodon, mutcodon, refaa, resnum, mutaa)
#   nucpos and resnum are all 1-based

def calc_mutation(genomepos,base,gene,genome):
  Rv,name = gene[2],gene[3]
  start,end = gene[0],gene[1] # annotation data is 0-based (once read in)
  seq = genome[start:end+1]
  n,x = len(seq),genomepos-start # x is pos of mut in gene (in bp)
  mut = base
  if gene[4]=='-':
    seq = reverse_complement(seq)
    x = n-x-1
    mut = complement.get(mut,'N')
  y = 3*int(x/3) # y is start of codon in gene (in bp)
  refcodon = seq[y:y+3]
  mutcodon = list(refcodon)
  mutcodon[x%3] = mut
  mutcodon = ''.join(mutcodon)
  c1 = codon.get(refcodon,"?")
  c2 = codon.get(mutcodon,"?")
  return (Rv,name,x+1,n,x%3+1,refcodon,mutcodon,c1,int(x/3)+1,c2)

complement = {'A':'T','T':'A','C':'G','G':'C'}

def reverse_complement(seq):
  s = list(seq)
  s.reverse()
  for i in range(len(s)):
    s[i] = complement.get(s[i],s[i]) # if unknown, leave as it, e.g > or !
  s = ''.join(s)
  return s

###############################

if len(sys.argv)<4: 
  sys.stderr.write("usage: python annotate_vcf.py <.vcf> <ref.fna> <ref.prot_table>\n")
  sys.exit(0)

genome = read_genome(sys.argv[2])
genes = read_genes(sys.argv[3])
hash = hash_genes(genes,genome)

# note: ignores Amb and lowCov sites (assume reference)

for line in open(sys.argv[1]):
  if line[0]=='#': print(line.rstrip()); continue
  w = line.rstrip().split('\t')
  pos = int(w[1])-1
  ref,nuc = w[3],w[4]
  gene = lookup_gene(pos,hash)
  orf = "intergenic" if gene==None else "%s/%s" % (gene[2],gene[3])
  temp = [orf]
  if w[6]=="PASS" and len(w[3])==1 and len(w[4])==1 and w[4]!='.': 
    temp += ['*','snp']
    if gene==None: temp.append("%s:%s>%s" % (pos+1,ref,nuc))
    else:
      result = calc_mutation(pos,nuc,gene,genome)
      temp.append("%s/%s:%s%s%s" % (result[0],result[1],result[7],result[8],result[9]))
  if w[6]=="PASS" and (len(w[3])>1 or len(w[4])>1): 
    temp.append('*')
    dn = len(w[4])-len(w[3])
    if dn>0:   temp += ['ins',"%s:+%s" % (orf,dn)]
    elif dn<0: temp += ['mdel',"%s:%s" % (orf,dn)] # consolidated del line
    else:      temp += ['subst',"%s:subst%s" % (orf,len(w[3]))]
  if "Del" in w[6]: temp += ['x','del'] # individual del sites are needed for "exclusion" in snp_dists.py

  cov,het = -1,-1
  info = w[7]
  x,y = info.find("DP="),info.find("BC=")
  if x!=-1 and y!=-1:
    DP = info[x+3:info.find(";",x)]
    BC = info[y+3:info.find(";",y)]
    cov = int(DP)
    nucs = BC.split(',')
    nucs = [int(x) for x in nucs]
    a,b = max(nucs),sum(nucs)
    het = 0 if b==0 else (b-a)/float(b)

  vals = w+[cov,round(het,3)]+temp
  print('\t'.join([str(x) for x in vals]))
  
