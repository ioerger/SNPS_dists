import sys

# based on table of SNPs in Coll et al (2014). https://www.nature.com/articles/ncomms5812

if len(sys.argv)<3:
  print("usage: python mtb_lineages_vcf.py <.vcf> coll_2014_LS_SNPs.txt")
  sys.exit(0)

SNPs = {}
for line in open(sys.argv[2]): # "coll_2014_LS_SNPs.txt"
  w = line.split()
  SNPs[w[1]] = w

lineages = []
for line in open(sys.argv[1]): # vcf file
  if line.startswith('#'): continue
  w = line.rstrip().split('\t')
  pos,kind = w[1],w[6]
  if pos in SNPs and kind=='PASS':
    v = SNPs[pos]
    mut = "%s/%s" % (w[3],w[4])
    lineage = v[0]
    if ("**" not in lineage and mut==SNPs[pos][3]) or ("**" in lineage and mut!=SNPs[pos][3]): print("# %s" % ('\t'.join(v))); lineages.append(lineage)

# the longest lineage name should be the most specific
print(sorted(lineages,reverse=True)[0])
