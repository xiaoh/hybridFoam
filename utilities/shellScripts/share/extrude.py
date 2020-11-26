#!/usr/bin/env python

# CONFIGURATION
infile = "blk.tmp"
outfile = "blk_extrude.tmp"
nz = 31
lz = 4

import numpy as np

A = np.loadtxt(infile, comments="#")
fout = open(outfile, 'wt')

l = 0
while True:
  n = A[l,:].prod()
  B = A[l+1:l+1+n,:]
  B.shape = B.shape[0], B.shape[1], 1
  B = B.repeat(nz, 2)
  B[:,2,:] = np.linspace(0,lz,nz).reshape((1,nz)).repeat(B.shape[0], 0)
  B = B.transpose((0,2,1)).reshape((-1,3))
  fout.write("%d %d %d\n"%(A[l,0], A[l,1], nz))
  np.savetxt(fout, B)
  l += n+1
  if l >= A.shape[0]:
    break
