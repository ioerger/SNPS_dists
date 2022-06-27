import sys

if len(sys.argv)<3: 
  print("usage: python snp_dists.py <sites_file> <treeorder>")
  sys.exit(0)

strains = []
for line in open(sys.argv[2]): # treeorder
  name = line.split()[0]
  #if "_" in name: name = name[:name.find("_")] # was for removing "_aby"; beware of P1_G11
  strains.append(name)

names,seqs = [],[]
#skip = 1
for line in open(sys.argv[1]): # sites
  #if skip>0: skip -= 1; continue
  if len(line)<20: continue # skip header line(s)
  w = line.split()
  name = w[0]
  #if "_" in name: name = name[:name.find("_")]
  names.append(name)
  seqs.append(w[1])

for name in strains:
  if name not in names: 
    sys.stderr.write("error: %s is in strains (treeorder) but not in names (sites)\n" % name)
    sys.exit(0)

for n,a in enumerate(strains):
  s = ""
  sys.stderr.write("%s: %s\n" % (n+1,a))
  s += "%10s " % a
  for b in strains:
    i,j = names.index(a),names.index(b)
    snps = 0
    for k in range(len(seqs[0])):
      if seqs[i][k]!=seqs[j][k]: snps += 1
    s += "%3d " % snps
  print(s)
  sys.stdout.flush()

