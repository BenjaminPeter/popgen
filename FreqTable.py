#!/usr/bin/env python
"""
    a FreqTable is a table where each row is a SNP and each
    column is a population
"""

import numpy as np
from Population import *


class FreqTable:
    def __init__(self, n_snp=0, sample_sizes=[], pops =[]):
        self.n_samples = len(sample_sizes)
        self.sample_sizes = sample_sizes
        self._data = np.zeros((n_snp, self.n_samples))
        self.pops = pops

    def __getitem__(self,x):
        return self._data.__getitem__(x)

    def __setitem__(self,*args):
        return self._data.__setitem__(*args)

    def __len__(self):
        return self._data.shape[0]

    @staticmethod
    def load(pop_file,data_file):
        self = FreqTable()

        self.pops = []
        for fh_pop_line in open(pop_file,"r"):
            self.pops.append( Population( line=fh_pop_line ) )
        for p in self.pops:
            self.sample_sizes.append( p.sample_size )
        self._data = np.loadtxt(data_file)
        self.n_snp, self.n_samples = self._data.shape
        return self

    #------- range expansion stuff ----------------

    def psi(self):
        n = self.n_samples
        psi = np.zeros( n * (n-1) / 2 )
        shared_snp = np.zeros( n * (n-1) / 2 )
        psiDict = {}
        col = 0
        for i in xrange(n-1):
            for j in xrange(i+1,n):
                psiDict[i,j] = col
                for snp in self:
                    if snp[i] and snp[j]:
                        shared_snp[ col ] += 1
                        psi[ col ] += snp[i] - snp[j]

                if shared_snp[ col ] == 0:
                    psi[col] = 0
                else:
                    psi [col] /= shared_snp[ col ]
                col +=1
        return psi, shared_snp



