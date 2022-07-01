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

sizes = []
for size,i in pairs:
  if size<1000: break
  print(h[i])
  print_seq(s[i])
  sizes += [size]*size

temp = list(filter(lambda x: x[0]>1000,pairs))
sys.stderr.write("contigs > 1kb:  %s\n" % len(temp))
sys.stderr.write("longest contig: %s bp\n" % pairs[0][0])
sys.stderr.write("total length:   %s bp\n" % (sum([x[0] for x in pairs])))
sys.stderr.write("N50:            %s bp\n" % sizes[int(len(sizes)/2)])
