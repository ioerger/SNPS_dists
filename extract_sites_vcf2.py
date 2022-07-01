import sys

genome_size = 0

def read_vcf2(fname):
  global genome_size
  muts = {}
  for line in open(fname):
    if line[0]=='#': continue
    w = line.rstrip().split('\t')
    if '*' not in line: continue
    co,ref,mut = int(w[1]),w[3],w[4]
    genome_size = max(co,genome_size) # max
    kind = w[11] # snp, ins, mdel (consolidated), or del (indiv nucs)
    info = "" if len(w)<13 else w[12]
    muts[co] = (kind,ref,mut,info)
  return muts

isolates = []
for line in open(sys.argv[1]): # file with vcf2 filenames
  w = line.rstrip().split()
  isolates.append(w[0])
vcf2_files = ["%s.vcf2" % x for x in isolates]

datasets = [] # hash tables of mutations
for fname in vcf2_files:
  datasets.append(read_vcf2(fname))
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
  vals = [vcf2_files[i],snp,ins,dels]
sites = sorted(list(sites))

# how many candidates sites are there? (not just snps)
#   exclude sites where there is an ins/del in any strain
#   exclude PPE/PGRS

snp_sites,bad_sites = [],[]
for site in sites:
  good = True
  for muts in datasets:
    if site in muts:
      if muts[site][0]!='snp': good = False
      if "PPE" in muts[site][-1]: good = False
      if "PGRS" in muts[site][-1]: good = False
  if good==True: snp_sites.append(site)
  else: bad_sites.append(site)
filtered = len(bad_sites)
candidates = genome_size-filtered
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
  for hash in datasets:
    if co not in hash: nucs.append("."); continue
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

