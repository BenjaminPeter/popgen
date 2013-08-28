import numpy as np
import unittest
from Haplotypes import Haplotypes
from Haplotypes import BijectiveDict

class TestMsmsData(unittest.TestCase):
    """
        my first try at an unittest. Loads a simple msms data set and
        runs some standard tests
    """
    def setUp(self):
        self.h = Haplotypes()
        self.h.read_ms_file("test/test.msms")

    def test_bijective_dict(self):
        d = BijectiveDict()
        for i in range(20):
            d[i] = i+30
        for i in range(20):
            self.assertEqual(d[i], i+30)
            self.assertEqual(d.r[i+30], i)
        l = len(d)
        del d[4]
        self.assertEqual(len(d), l-1)
        self.assertEqual(len(d.r), l-1)

    def test_subsetters(self):
        """
            test if the selectors work properly
        """
        subset = self.h[(10,20),(3,5)]
        self.assertTrue(subset.shape == (2,10))

        self.h._default_coordinate_system = 'pos'
        subset = self.h[(.2,.6),(5,8)]


    def test_rec(self):
        self.h
        pass
    def test_sfs_freqt(self):
        """
            test the creation of the SFS and freqtable
        """
        return
        sfs,freqt = self.h.make_SFS()
        self.assertEqual(len(freqt),sum(sfs.sfs))
        for i,snp in enumerate(self.h):
            self.assertEqual(freqt[i],sum(snp))

        for i, freq in enumerate(sfs.sfs):
            self.assertEqual(freq, sum(freqt._data == i))

        self.h.default_coordinate_system = "id"
        print self.h.default_coordinate_system

if __name__ == '__main__':
    unittest.main()
