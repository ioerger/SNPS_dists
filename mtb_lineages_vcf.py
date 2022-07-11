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
problematic = False
for line in open(sys.argv[1]): # vcf file
  if line.startswith('#'): continue
  w = line.rstrip().split('\t')
  pos,kind = w[1],w[6]
  if pos in SNPs:
    v = SNPs[pos]
    mut = "%s/%s" % (w[3],w[4])
    lineage = v[0]
    if lineage=="lineage4.5" and "lineage1" in lineages: print("# skipping analysis of lineage4.5 in presence of SNP for lineage1"); continue
    if kind=='PASS':
      if ("**" not in lineage and mut==SNPs[pos][3]) or ("**" in lineage and mut!=SNPs[pos][3]):
         print("# %s" % ('\t'.join(v))); lineages.append(lineage)
    elif "LowCov" in kind:
      print("# %s" % ('\t'.join(['warning:',lineage,pos,mut])))
      print("#   vcf: %s" % (line[:-1]))
      print("#   Low-coverage site - could be involved in an indel")
      problematic = True
    elif "Amb" in kind:
      print("# %s" % ('\t'.join(['warning:',lineage,pos,mut])))
      print("#   vcf: %s" % (line[:-1]))
      print("#   Ambiguous lineage-specific sites could indicate a mixed isolate")
      problematic = True

if problematic: print "unclear" # this is different from 'unrecognized'
elif len(lineages)==0: print("unrecognized")
else:
  baselin = lineages[0][7] # digit, like 1-7 or B (for BOV)
  allsame = True
  for lin in lineages:
    if lin[7]!=baselin: allsame = False
  if allsame==False: print("inconsistent") # SNPs matching multiple lineages
  else: print(sorted(lineages,reverse=True)[0]) # the longest lineage name should be the most specific
