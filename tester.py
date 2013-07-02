import numpy as np
from Haplotypes import Haplotypes
from types import GeneratorType

h = Haplotypes()
testPos=(46470000,46480000)
testHaps=np.arange(40)
h.read_mbs_file("../../selectiveSweep/data/originalData/tib_v3.mbs")
idGen, haps = h._get_default_location(pos=testPos)

#get some data subset
data = h.get_SNP_data(pos=testPos,haplotypes = testHaps)
