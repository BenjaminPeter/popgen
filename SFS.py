#!/usr/bin/env python

import numpy as np
class SFS2D():
    def __init__(self,n):
        self.sfs = np.zeros((n[0]+1,n[1]+1), dtype=int)
        self.n_hap = n

        self.sfs1d = None,None


    @staticmethod
    def loadFromHaplotypes(h, p0, p1, **kwargs):
        """
        """
        id, haplotypes = h._get_default_location(**kwargs)
        if type(p0[0]) == bool or type(p0[0]) is np.bool_:
            n0 = sum(p0)
        else:
            n0 = len(p0)
        if type(p1[0]) == bool or type(p1[0]) is np.bool_:
            n1 = sum(p1)
        else:
            n1 = len(p1)

        sfs = SFS2D((n0,n1))
        for i in id:
            snp = h.get_SNP_data(id=i, haplotypes=haplotypes)
            p = sfs.addSNP(snp,p0,p1)

        return sfs

    def __getitem__(self,x):
        return self.sfs.__getitem__(x)

    
    def make_1d_sfs(self):
        self.sfs1d = np.sum(self.sfs,0), np.sum(self.sfs,1)
        return self.sfs1d

    def addSNP(self,snp,pop1,pop2):
        """snp is the genotype for each haplotype
        pop1,pop2 are arrays/lists s.t. snp[pop1] 
        is the number of alleles in pop1 """
        p1 = sum( snp[pop1] )
        p2 = sum( snp[pop2] )
        self.sfs[p1,p2] += 1
        return p1,p2

    def Pi(self):
        if self.sfs1d is None:
            self.make_1d_sfs()
        return self.sfs1d[0].Pi(), self.sfs1d[1].Pi()

    def S(self):
        if self.sfs1d is None:
            self.make_1d_sfs()
        return self.sfs1d[0].S(), self.sfs1d[1].S()

    def H(self):
        if self.sfs1d is None:
            self.make_1d_sfs()
        return self.sfs1d[0].H(), self.sfs1d[1].H()

    def wattersons_theta(self):
        if self.sfs1d is None:
            self.make_1d_sfs()
        return self.sfs1d[0].wattersons_theta(),\
                self.sfs1d[1].wattersons_theta()

    def theta_h(self):
        if self.sfs1d is None:
            self.make_1d_sfs()
        return self.sfs1d[0].theta_h(), self.sfs1d[1].theta_h()

    def fay_wu_h(self):
        if self.sfs1d is None:
            self.make_1d_sfs()
        return self.sfs1d[0].fay_wu_h(), self.sfs1d[1].fay_wu_h()

    def tajimas_d(self):
        if self.sfs1d is None:
            self.make_1d_sfs()
        return self.sfs1d[0].tajimas_d(), self.sfs1d[1].tajimas_d()

    def __getitem__(self,x):
        return self.sfs.__getitem__(x)



class SFS():
    def __init__(self, n):
        """Initializes an empty 1D SFS with a sample size of n
            verbose: should debug stuff be printed
        """
        self.sfs = np.zeros(n+1,dtype=int)
        self.n_hap = n


        self._pi = None
        self._td = None
        self._S = None
        self._FWH = None
        self._ns = None
        self._theta_h = None

    def __getitem__(self,pos):
        return self.sfs[pos]

    def __setitem__(self,pos, data):
        self.sfs[pos] = data
    def addSNP(self,snp):
        p = sum( snp )
        self.sfs[p] += 1
        return p
    
    @staticmethod
    def loadFromHaplotypes(h, **kwargs):
        """creates a new SFS from a haplotype object"""
        raise ValueError("DONT USE")
        id, haplotypes = h._get_default_location(**kwargs)
        nHap = len(haplotypes)

        sfs=SFS(nHap)

        """freqs here should become a freqTable object"""
        freqs=[]

        for i in id:
            snp = h.get_SNP_data(id=i,haplotypes=haplotypes)
            p = sfs.addSNP(snp)
            #freqs.append( float(p) )
        return sfs, freqs


    def get_H(self):
        """calculates heterozygosity"""
        return self.Pi() / self.S()
    get_Heterozygosity = get_H

    def Pi(self):
        if self._pi is not None:
            return self._pi

        pi = 0.0
        n = self.n_hap
        nnm1 = n*(n-1.)/2.
        for i,site in enumerate(self.sfs):
            pi += site * i * (n-i)

        self._pi = pi
        return pi / nnm1
    n_pw_differences = Pi
    tajimas_estimator = Pi

    def S(self):
        if self._S is not None:
            return self._S

        S = sum(self.sfs[1:self.n_hap])
        self._S = S
        return S
    n_segregating_sites = S

    def wattersons_theta(self):
        hs = 0.0
        for i in xrange(1,self.n_hap):
            hs += 1./i

        return self.S() / hs
    wattersons_estimator = wattersons_theta

    def theta_h(self):
        """calculates fay and wu's theta H"""
        if self._theta_h is not None:
            return self._S
        theta_h = 0.0
        n = self.n_hap
        nnm1 = n*(n-1.)/2.
        for i,site in enumerate(self.sfs):
            pi += site * i * i

    def fay_wu_h(self):
        """add some checks about polarized data"""
        if self._FWH is not None:
            return self._FWH
        self._FWH = self.Pi() - self.theta_h()
        return self._FWH

    def tajimas_d(self):
        if self._td is not None:
            return self._td

        a1=0.0
        n= self.n_hap
        for i in range(1,n):
            a1 += 1./i
        #first, calculate coefficients a2,e1 and e2
        a2 = 0.
        for i in range(1,n):
            a2 += 1. / ( i * i )
        e1 = 1./a1 * ((n+1.)/(3.*(n-1.))-1/a1)
        e2 = 1./(a1*a1+a2) * (((2.*(n*n+n-3.))/(9.*n*(n-1))) - ((n+2.)/(n*a1))+a2/(a1*a1))
        vartD=e1*S+e2*S*(S-1.)

        self._td = ( self.Pi() - self.wattersons_theta() )/np.sqrt(vartD)

    def __repr__(self):
        return self.sfs.__repr__()
