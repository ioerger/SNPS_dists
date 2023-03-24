import sys

genome_size = 0

def read_vcf2(fname):
  global genome_size
  muts,nonpassing = {},{}
  for line in open(fname):
    if line[0]=='#': continue
    w = line.rstrip().split('\t')
    if w[6]!="PASS": nonpassing[int(w[1])] = 1
    if len(w)<12 or w[11] not in '*x': continue
    co,ref,mut = int(w[1]),w[3],w[4]
    genome_size = max(co,genome_size) # max
    kind = w[12] # snp, ins, mdel (consolidated), or idel (indiv nucs)
    info = "" if len(w)<14 else w[13]
    muts[co] = (kind,ref,mut,info)
  return muts,nonpassing

if len(sys.argv)<3: 
  print("usage: python extract_sites_vcf2.py <file_of_strain_names> <prot_table> [-vertical]")
  sys.exit(0)

isolates = []
for line in open(sys.argv[1]): # file with vcf2 filenames
  w = line.rstrip().split()
  isolates.append(w[0])
vcf2_files = ["%s.vcf2" % x for x in isolates]

datasets = [] # hash tables of mutations
nonpassingsites = []
for fname in vcf2_files:
  muts,nonpassing = read_vcf2(fname)
  nonpassingsites.append(nonpassing)
  datasets.append(muts)
  sys.stderr.write("done reading %s\n" % fname)
  sys.stderr.flush()

# make a list of all sites in all strains
# filter out any site that is an indel in any strain

sites = set()
for i,muts in enumerate(datasets): 
  sites = sites.union(muts.keys())
  n = len(muts) # snps and indels; no good - this counts dels as individual nucs
  snp = len(list(filter(lambda x: x[0]=='snp',muts.values()))) # just snps
  ins = len(list(filter(lambda x: x[0]=='ins',muts.values()))) 
  dels = len(list(filter(lambda x: x[0]=='mdel',muts.values()))) # consolidated dels
  vals = [vcf2_files[i],snp,ins,dels] # is this used anymore?
sites = sorted(list(sites))

# how many candidates sites are there? (not just snps)
#   exclude sites where there is an ins/del in any strain
#   exclude PPE/PGRS

EXCLUDE_PPE = True
EXCLUDE_PGRS = True
EXCLUDE_INDELS = True
if "+PPE" in sys.argv: EXCLUDE_PPE = False
if "+indels" in sys.argv: EXCLUDE_INDELS = False

snp_sites,bad_sites = [],[]
for site in sites:
  good = True
  for muts in datasets:
    if site in muts:
      if EXCLUDE_INDELS==False and muts[site][0] not in ['snp','ins','mdel']: good = False
      elif EXCLUDE_INDELS==True and muts[site][0]!='snp': good = False
      if EXCLUDE_PPE==True and "PPE" in muts[site][-1]: good = False
      if EXCLUDE_PGRS==True and "PGRS" in muts[site][-1]: good = False
  if good==True: snp_sites.append(site)
  else: bad_sites.append(site)
filtered = len(bad_sites) # includes indels, but only among sites with a mutation

bad_sites_hash = {}
for i in bad_sites: bad_sites_hash[i] = 1
for line in open(sys.argv[2]):
  w = line.rstrip().split('\t')
  if (EXCLUDE_PPE and "PPE" in w[7]) or (EXCLUDE_PGRS and "PGRS" in w[7]): 
    for j in range(int(w[1]),int(w[2])+1): bad_sites_hash[j] = 1
candidates = genome_size-len(bad_sites_hash.keys())
rate = len(snp_sites)/float(candidates)

#sys.stderr.write("num of isolates: %s, ref genome size: %s, sites filtered out: %s, candidate sites: %s, coverage: %0.1f%%, SNP sites: %s, \n" % (len(datasets),genome_size,filtered,candidates,candidates*100/float(genome_size),len(snp_sites)))

sys.stderr.write("number of isolates:  %s\n" % len(datasets))
sys.stderr.write("ref genome size:     %s\n" % genome_size)
sys.stderr.write("sites filtered out:  %s\n" % filtered)
sys.stderr.write("candidate sites:     %s\n" % candidates)
sys.stderr.write("coverage:            %0.2f%%\n" % (candidates*100/float(genome_size)))
sys.stderr.write("allelic (SNP) sites: %s\n" % (len(snp_sites)))
sys.stderr.write("allelic frequency:   %0.6f (%0.2f SNP sites per 1kb)\n" % (rate,rate*1000))

snplists = []
for i,hash in enumerate(datasets):
  subset = []
  for coord,mut in hash.items():
    if coord in snp_sites: subset.append((coord,mut))
  snplists.append(subset)

VERTICAL = False
if "-vertical" in sys.argv: VERTICAL = True

snp_sites.sort()
n = len(datasets)
allnucs = []
if VERTICAL==True: print('\t'.join("coord H37Rv".split()+isolates))
for i in range(len(snp_sites)):
  co = snp_sites[i]
  ref,info,nucs = None,None,[]
  #for hash in datasets:
  #  if co not in hash: nucs.append("."); continue
  for hash,nonpassing in zip(datasets,nonpassingsites):
    if co not in hash:
      if co in nonpassing: nucs.append("?")
      else: nucs.append("."); 
      continue
    # rewrite ref and mut with each occurence; assume translation is the same for all
    kind,ref,mut,info = hash[co] 
    nucs.append(hash[co][2])
  vals = [co,ref]+nucs+[info]
  if VERTICAL==True: print('\t'.join([str(x) for x in vals])) # show table of sites as rows)
  nucs = [ref if x=="." else x for x in nucs]
  allnucs.append(nucs)
if VERTICAL==True: sys.exit(0)

a,b = len(allnucs),len(allnucs[0])
#rotated = [""]*b
#for i in range(a):
#  for j in range(b): rotated[j] += allnucs[i][j]
rotated = ["".join([x[i] for x in allnucs for i in range(b)])]
print("xread")
print("%s %s" % (a,b))
for i in range(b):
  #print("> %s" % vcf2_files[i])
  nucs = "".join([row[i] for row in allnucs])
  #print("%-10s %s" % (vcf2_files[i][:10],nucs))
  print("%s %s" % (isolates[i],nucs))

