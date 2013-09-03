import numpy as np
import unittest
from Haplotypes import *

class TestMsmsData(unittest.TestCase):
    """
        my first try at an unittest. Loads a simple msms data set and
        runs some standard tests
    """
    def setUp(self):
        self.h = Haplotypes()
        self.h.read_ms_file("test/test.msms")

    def test_coords(self):
        d = Coords(np.arange(0,1,.01))
        self.assertTrue(np.all(d.f == sorted(d.f)))
        for k,v in d.r.d.items():
            self.assertEqual(k,d[v])

        #setting
        d[4] = 0.84
        self.assertEqual(d.f, sorted(d.f))
    
        #inserting
        d.insert(34,.345)
        self.assertTrue(np.all(d.f == sorted(d.f)))
        for k,v in d.r.d.items():
            self.assertEqual(k,d[v])

        #append
        d.append(1.01)
        self.assertTrue(np.all(d.f == sorted(d.f)))
        for k,v in d.r.d.items():
            self.assertEqual(k,d[v])

    def test_subsetters(self):
        """
            test if the selectors work properly
        """
        subset = self.h[10:20,3:7]._data
        self.assertTrue(subset.shape == (4,10))

        subset = self.h[10:20]._data
        self.assertTrue(subset.shape == (self.h.n_hap,10))
        
        subset = self.h[:,4:19]._data
        self.assertTrue(subset.shape == (15,self.h.n_snp))

        subset = self.h[:,[4,5,13]]._data
        self.assertTrue(subset.shape == (3,self.h.n_snp))

        subset = self.h[:,23]._data
        self.assertTrue(subset.shape == (1,37))

        subset = self.h[10]._data
        self.assertTrue(subset.shape == (self.h.n_hap,1))

        subset = self.h[10,:]._data
        self.assertTrue(subset.shape == (self.h.n_hap,1))

        subset = self.h[10,:10]._data
        self.assertTrue(subset.shape == (10,1))

    def test_subsetters_pos(self):
        """ tests if subsetters work using genetic positions"""
        self.h.default_coordinate_system='pos'

        subset = self.h[:,[23,43]]._data
        self.assertTrue(subset.shape == (2,37))

        subset = self.h[:.3,43]._data
        self.assertEqual(subset.shape[1],sum(self.h.coords['pos'].f<0.3))

        subset = self.h[.2:,23]._data
        self.assertEqual(subset.shape[1],sum(self.h.coords['pos'].f>=0.2))

        subset = self.h[0.01809]._data
        self.assertEqual(subset.shape,(100,1))

        subset = self.h[[0.04667,0.70350]]._data
        self.assertEqual(subset.shape,(100,2))

    def test_sfs_freqt(self):
        """
            test the creation of the SFS and freqtable
        """
        sfs,freqt = self.h.make_SFS()
        self.assertEqual(len(freqt),sum(sfs.sfs))
        for i,snp in enumerate(self.h._data.transpose()):
            self.assertEqual(freqt[i],sum(snp))

        for i, freq in enumerate(sfs.sfs):
            self.assertEqual(freq, sum(freqt._data == i))

        self.h.default_coordinate_system = "id"
        print self.h.default_coordinate_system

    def test_EHH_IHS(self):
        self.assertEqual(len(self.h._unique_haplotypes()),21)
        self.assertNotEqual(self.h.EHH((3,6),id=11,coords='id'),self.h.EHH((.3,.6),id=11))


if __name__ == '__main__':
    unittest.main()
