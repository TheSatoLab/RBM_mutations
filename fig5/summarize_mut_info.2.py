#!/usr/bin/env python

import sys, re
argvs = sys.argv

f = open(argvs[1])

f.__next__()

##2 columns
print("Id\tmut")
for line in f:
  line = line.strip().split("\t")
  Id = line[4] #5th col
  mut_line = line[16].replace('(','').replace(')','') #17th col
  mut_l = mut_line.split(',')
  for mut in mut_l:
    print("\t".join([Id,mut]))

##3columns
# print("gene\tmut\tId")
# for line in f:
#   line = line.strip().split("\t")
#   Id = line[2]
#   mut_line = line[14].replace('(','').replace(')','')
#   if len(mut_line) > 0:
#     mut_l = mut_line.split(',')
#     for mut in mut_l:
#       mut_p = mut.split("_")
#       gene = mut_p[0]
#       mut_d = mut_p[1]
#       print("\t".join([gene,mut_d,Id]))
