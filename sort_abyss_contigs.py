import sys

def read_fasta(filename):
  headers,seqs = [],[]
  seq = ""
  for line in open(filename):
    line = line.rstrip()
    if len(line)==0: continue
    if line[0]==">":
      headers.append(line)
      if seq!="": seqs.append(seq)
      seq = ""
    else: seq += line
  seqs.append(seq)
  return headers,seqs

def print_seq(s,W=70):
  i,n = 0,len(s)
  while i<n:
    print(s[i:i+W])
    i += W

##################

h,s = read_fasta(sys.argv[1])

pairs = []
for i in range(len(h)):
  size = int(h[i].split()[1]) # length of contig
  pairs.append((size,i))

pairs.sort(reverse=True)

for size,i in pairs:
  if size<1000: break
  print(h[i])
  print_seq(s[i])

temp = list(filter(lambda x: x[0]>1000,pairs))
sys.stderr.write("%s contigs of size > 1kb\n" % len(temp))
sys.stderr.write("longest contig: %sbp\n" % pairs[0][0])
sys.stderr.write("total length: %sbp\n" % (sum([x[0] for x in pairs])))
