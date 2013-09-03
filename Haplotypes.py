#!/usr/bin/env python


import numpy as np
from types import GeneratorType
from SFS import SFS
from FreqTable import FreqTable

class Coords(object):
    """
        handles stuff like id <==> position relation
    """
    class RevDict(object):
        def __init__(self, parent):
            self.d = dict()
            self.parent = parent
        def __getitem__(self,x):
            if isinstance(x,slice):
                xstart, xstop = x.start,x.stop
                if x.step is not None:
                    raise ValueError("Slice not supported")
                if x.start is None:
                    xstart = 0
                if x.stop is None:
                    xstop = np.finfo('d').max 

                id0 = self.get_closest_SNP_from_pos(xstart-1e-12,1)
                id1 = self.get_closest_SNP_from_pos(xstop-1e-12,-1)
                return np.arange(id0,id1+1)
            elif isinstance(x,tuple):
                id0 = self.get_closest_SNP_from_pos(x[0],1)
                id1 = self.get_closest_SNP_from_pos(x[1],-1)
                return np.arange(id0,id1)
            elif isinstance(x,list) or isinstance(x,np.ndarray):
                results = [self.get_closest_SNP_from_pos(i) for i in x ]
                return np.array(results)
            else:
                id = self.get_closest_SNP_from_pos(x)
                return np.array([id])
        def __setitem__(self,x,y):
            self.d.__setitem__(x,y)

        def get_closest_SNP_from_pos(self,pos,direction=0):
            """gives a SNP close to the position. If direction==None or direction
            ==0, the closest SNP
            is returned. is direction==1, the next larger SNP is returned, if
            direction== =1, the next SNP with lower position is returned """

            try:
                s = self.d[pos]
                return s
            except:
                pass

            #do binary search: O(lg n)
            segs= self.parent.f
            l = len(segs)
            while l>1:
                pivot = l/2
                if pos> segs[pivot]: 
                    segs = segs[pivot:]
                else:
                    segs = segs[:pivot]
                l= len(segs)
            lower = segs[0]
            lPos = self.d[lower]
            if lPos ==0:
                if direction==-1:
                    raise ValueError('no smaller SNP than this')
                elif pos < self.parent.f[0]:
                    return 0 
                else:
                    return 1


            uPos = lPos+1
            if uPos ==len(self.parent.f):
                if direction==1:
                    raise ValueError('no SNP genotyped at larger pos')
                else:
                    return lPos
            upper = self.parent.f[uPos]

            #print lPos, uPos, "diff: ", pos-lower, upper-pos
            if abs(pos - lower) > abs(upper-pos) and direction != -1:
                return uPos
            elif direction != 1:
                return lPos
            return uPos

        def clear(self):
            self.d.clear()

    def __init__(self, *args, **kwargs):
        self.f = np.array(list(*args, **kwargs))
        self.r = self.RevDict(self)
        self.update_r_dict()

    def update_r_dict(self):
        """
            updates the dict, done when an item is changed or inserted
            this is of course rather inefficient and O(n)
        """
        self.r.clear()
        for i,e in enumerate(self.f):
            self.r[e] = i

    def __setitem__(self,x,y):
        """
            the most general setter. As the list is guaranteed to be
            ordered, the list will be sorted and the dict updated
        """
        self.f.__setitem__(x,y)
        self.r.__setitem__(y,x)
        self.f = sorted(self.f)
        self.update_r_dict()

    def __getitem__(self,x):
        return self.f.__getitem__(x)

    def __len__(self):
        return self.f.__len__()

    def __str__(self):
        return self.f.__str__()

    def __repr__(self):
        return self.f.__repr__()

    def insert(self,x,y):
        self.f.insert(x,y)
        self.r[y] = x
        if x < len(self)-1:
            assert self.f[x] < self.f[x+1]
        if x > 0:
            assert self.f[x-1] < self.f[x]
        self.update_r_dict()

    def append(self,x):
        return self.insert(len(self),x)

class IDCoords(Coords):
    """
        a simple "identity subclass that handles id retrival
    """
    def __init__(self, *args, **kwargs):
        self.r = self
        return

    def __getitem__(self, x):
        return x

    def __setitem__(self, *args):
        return

    def update_r_dict(self):
        return 
    
    def __repr__(self):
        return "obj"



class Haplotypes(object):
    """
        this class characterizes a set of haplotypes 
        it is intended for small to intermediate data sets,
        such that all the haplotypes can be kept in memory
        simultaneously. This is the data structure that saves most of the information, vs the SFS, where the linkage information is lost and the FreqTable, where we care about many populations

        As a data structure, I'll use a np.array to store all the data,
        For the genomic positions of the SNP, I'll use a 1D masked array
        that has all SNP positions and a dict[genomic position in bp] -> snp id
        for the recombination position of the SNP, I'll use a 1D masked array that has all SNP positions as a float
        pos or position refers to the genomic position, while id refers
        to the position in the array

        data will encoeded MS style, i.e 0 for ancestral, 1,2,3 for 
        derived alleles (support for 2,3 nyi). Bases i.e. ACGT will be
        kept in a separate numpy array

        haplotypes is a dict[individualID] = rows in array
        other attributes:
            polarized:      is the ancestral state known?
    """
    def __init__(self,file=None,verbose=False):
        """if file is given, it loads the file assuming mbs data"""
        """if verbose, some stuff is printed out"""
        """ should be extended to read files guessed by ending"""
        if file != None:
            self.read_mbs_file(file)

        self.reset()

        self.verbose = verbose
        
    def reset(self):
        """
            this function clears all the data that is present in the
            object, intended for cleaning up loading files, etc.
        """
        #dict id <=> genomic position in bp
        self._seg_sites = Coords()
        #dict id <=> genomic position in rec distance
        self.coords = dict()
        self.coords['pos'] = self._seg_sites
        self.coords['bp']  = self._seg_sites
        self.coords['rec'] = self._seg_sites

        self.sfs = None
        self._data = None


        self.polarized = None
        self.alleles = None

        #dict hid <=> individual name
        self.individual_names = None

        self._default_coordinate_system = 'id'
        self._default_individual_selector = 'hid'

    def __len__(self):
        return self._data.shape[1]

#----------- Coordination and Indexing --------------------------
    def __getitem__(self,*args):
        """
            access to data should be as follows:
                - first coordinate is the SNP locations, using the 
                current coordinate system
                - second coordinate (if present) is the individuals
                e.g. if the coordinate system is pos, then h[3000:4000,10:20] should give all SNP betweehn position 3000 and 4000 from individuals with ids 10-19

        """
        if isinstance(args[0],tuple):
            id,haps = self._get_default(args[0][0], args[0][1])
        else:
            id,haps = self._get_default(args[0])
        if len(self._data[haps].shape) == 1:
            return self.subset(id,haps)
            return self._data[haps][id]
        return self.subset(id,haps)
        return self._data[haps][:,id]

    def _get_default(self, *args, **kwargs):
        """
            parses haplotype and location arguments. The first default
            arg in passed to _get_default_location, the second one is
            passed to _get_default_individuals. kwargs are passed to both. To ensure this works kwargs should be exclusive for the two.
        """
        if len(args)>0:
            loc = self._get_default_location(args[0], **kwargs)
        else:
            loc = self._get_default_location(**kwargs)
        if len(args)>1:
            haps = self._get_default_individuals(args[1], **kwargs)
        else:
            haps = self._get_default_individuals(**kwargs)

        return loc, haps

    def _get_default_individuals(self, *args, **kwargs):
        """
            function to handle kwargs for individuals. If none are 
            given, all individuals are used. Otherwise it takes a tuple
            of populations and takes those individuals, or a bunch of
            individual ids/names

            supported keywords are:
                pops: population objects on the current data set
                hid: the haplotype ids (numbered from 0 to n)
                name: the sample name of the samples, if defined
        """

        #the default here is to return all inds
        haplotypes = np.arange( self.n_hap )
        if len(args)>0:
            kwargs[self._default_individual_selector] = args[0]
        if len(args)>1:
            raise ValueError("Too many unnamed args")
        keys = kwargs.keys()
        if 'pops' in keys:
            haplotypes = []
            try:
                for pop in kwargs['pops']:
                    for sample in pop.samples:
                        haplotypes.append(sample)
            except:
                for sample in kwargs['pops']:
                    haplotypes.append(sample)


        elif 'hid' in keys:
            haplotypes=kwargs['hid']
        elif 'name' in keys:
            haplotypes=[self.individual_names.r[name] \
                        for name in kwargs['name']]


        return haplotypes

    def _get_default_location(self,*args, **kwargs):
        """
            function to handle kwargs for locations. If none are given,
            statistics are calculated over the entire data set and
            all individuals. otherwise, the following objects are
            supported:
                single number: calculate stuff for a single site
                tuple: beginning and end
                list or array: these sites
                slice: nyi

            supported keywords are:
                id: the id'th SNP
                pos: the genomic position is basepairs
                rec: the genomic position in recomb. units
        """

        cType = None
        if len(args)>0:
            kwargs[self.default_coordinate_system] = args[0]
            cType = self.default_coordinate_system
        if len(args)>1:
            raise ValueError("Too many unnamed args")
        keys = kwargs.keys()

        if sum(('id' in keys,'pos' in keys, 'rec' in keys))>1:
            raise ValueError("Error: multiples of id, pos and "+\
                             "rec specified")
        elif 'id' in keys:
            cType = 'id'
        else:
            for key in keys:
                if key in self.coords.keys():
                    cType = key
                    break

        #default if no coords are given
        if cType is None:
            if  'singleSite' in kwargs and kwargs['singleSite']:
                coords = [self.selId]
            else:
                coords = np.arange( len(self))
            cType = 'id'

        else:
            coords = kwargs[cType]


        if cType != 'id':
            coords = self.coords[cType].r[coords]


        return coords

    def subset(self, *args, **kwargs):
        """
            returns a new Haplotype object that displays a subset of the 
            original data. Rules whether data is copied are taken over from numpy, i.e. fancy indexing results in a copy, simple indexing doesn't
        """
        id, haps = self._get_default(*args, **kwargs)
        hnew = Haplotypes()
        hnew.reset()

        if self.polarized is not None:
            hnew.polarized = self.polarized[id]

        if self.alleles is not None:
            hnew.alleles = self.alleles[id]

        if self.individual_names is not None:
            hnew.individual_names = self.individual_names[haps]
        hnew.default_coordinate_system = self.default_coordinate_system
        hnew.default_individual_selector = self.default_individual_selector
        
        hnew._data = self._data[haps]
        if len(hnew._data.shape) == 1:
            hnew._data = hnew._data[id]
            hnew._data = hnew._data[np.newaxis]
        else:
            hnew._data = hnew._data[:,id]

        if len(hnew._data.shape) == 1:
            hnew._data = hnew._data[np.newaxis].transpose()

        for c in self.coords:
            if c == 'bp' or c == 'id':
                continue
            try:
                hnew.coords[c] = Coords(self.coords[c][id])
            except:
                hnew.coords[c] = Coords(self.coords[c][[id]])
        hnew.coords['bp'] = hnew.coords['pos']
        
        return hnew

#-----------  --------------------------
    def get_SNP_data(self, *args, **kwargs):
        """get the SNP data from either a position(default) or ID, returns the numpy array"""
        id, haps = self._get_default(singleSite=True, *args, **kwargs)
        return self._data[haps][:,id]

    def set_selected_site(self,*args, **kwargs):
        """
            function to assign a site to be the selected site.
            The main effect of this is for default behaviour of
            statistic when nothing else is declared
        """
        id = self._get_default_location(singleSite=True, *args, **kwargs)
        self.sel_id=id
        self.sel_coords = dict()
        for k in self.coords.keys():
            self.sel_coords[k] = self.coords[k][self.sel_id]
        #self.sel_site_data=self.get_SNP_data(id = self.sel_id)

    def get_haplotypes_with_derived_allel(self,*args, **kwargs):
        """gets all the haplotypes with the derived allel at the given snp. If id==None, the selected site is assumed"""
        data = self.get_SNP_data(**kwargs)
        return np.where(data==1)[0]

    def get_haplotypes_with_ancestral_allel(self,*args, **kwargs):
        """gets all the haplotypes with the ancestral allel at the given snp. If id==None, the selected site is assumed"""
        data = self.get_SNP_data(**kwargs)
        return np.where(data==0)[0]

    @property
    def default_coordinate_system(self):
        return self._default_coordinate_system

    @default_coordinate_system.setter
    def default_coordinate_system(self, value):
        if value not in ("id","rec","pos") and value not in self.coords.keys():
            raise ValueError("Tried to set coordinate"+\
                             " system to %s"%value)
        self._default_coordinate_system = value

    @property
    def default_individual_selector(self):
        return self._default_individual_selector

    @default_individual_selector.setter
    def default_individual_selector(self, value):
        if value not in ("hid","pop","name"):
            raise ValueError("Tried to set selector"+\
                             " system to %s"%value)
        self._default_individual_selector = value


    @property
    def n_hap(self):
        """returns the number of haplotypes in the data set"""
        return self._data.shape[0]

    @property
    def n_snp(self):
        """returns the number of SNP in the data sets"""
        return self._data.shape[1]

    @property
    def data(self):
        return self._data


#-----------------------------IHS/EHH----------------------------
    def _unique_haplotypes(self,*args,**kwargs):
        """for the data given in data, computes the unique haplotype, returns a dict"""
        id,haplotypes=self._get_default(*args, **kwargs)
        data=self[id,haplotypes]

        n_hap=data.n_hap
        unique=dict()
        #if we only have 1 segregating site, return allele frequencies
        if data.n_snp == 1:
            unique[0] = sum (data.data == 0)
            unique[1] = sum (data.data == 1)
            return unique

        #else, change to string... 
        for i in range(n_hap):
            s = tuple(j for j in data.data[i])
            if s in unique:
                unique[s]+=1
            else:
                unique[s] = 1

        return unique

    def _EHH_single(self, x, coords='bp', onesided=0, *args,**kwargs):
        """gets EHH for a single reference position x. Use EHH instead."""

        id , haplotypes=self._get_default( singleSite=True, *args,**kwargs)
        """call get EHH instead"""
        if haplotypes.shape == (0,):
            return 0 
        if onesided > 0:
            hapData=self.get_SNP_data(pos=(pos, pos+x), haplotypes=haplotypes)
        elif onesided < 0:
            hapData=self.get_SNP_data(pos=(pos-x,pos), haplotypes=haplotypes)
        else:
            pos = {coords:slice(self.coords[coords][id]-x, 
                                self.coords[coords][id]+x)}
            pos['hid'] = haplotypes
            
        return self._EHH_general(**pos)

    def EHH(self,x,onesided=0,*args, **kwargs):
        """calculates EHH in the range SelectedSite +/- x sites
                - onesided: if negative, only upstream allels are considered, if
                positive, only downstream allels are considered, otherwise both
                are considered(default)
                - id: the id of the SNP, if None, the selected SNP is assumed
                - haplotypes: if only a subset of haplotypes should be used:
                    default=all haplotypes
                    - verbose: adds some output"""
        ehh=[]
        if isinstance(x,tuple) or isinstance(x,np.ndarray):
            for xx in x:
                ehh.append(self._EHH_single(xx, onesided=onesided, 
                                            *args, **kwargs))
            return ehh
        else:
            return self._EHH_single(x,onesided=onesided, *args, **kwargs)

    def _EHH_general(self,*args, **kwargs):
        """don't call this. call EHH instead
            calculates EHH on the haplotypes given in hapData:
            hapData is expected to be an np.array of ints/bools where each row is a haplotype and each column is an individual
        """
        uniqueHaplotypes=self._unique_haplotypes(*args, **kwargs)
        cnt = uniqueHaplotypes.values()
        
        #get Denominator 
        denom = sum(cnt) * (sum(cnt)-1) 
        if denom == 0:
            return 0
            
        try:
            cnt = [cnt.count(i+1) for i in range(max(cnt))]
        except:
            raise ValueError(str(uniqueHaplotypes))


        nominator=0
        #print "found %i haplotypes at freq %i"% (cnt[0],1)
        for i in range(1,len(cnt)):
            nominator+=cnt[i]*i*(i+1)
            #print "found %i haplotypes at freq %i"% (cnt[i],i)

        return float(nominator) / denom
    
    def _integrate_EHH(self,dir,threshold,maxDist=None,
                         interpolation=False, verbose=False, **kwargs):
        """calculates the integrated EHH
            - id = id of the SNP to start at
            - dir = direction, if negative, upstream, if positive, downstream
            - ind = the haplotypes to consider
            - threshold = threshold where to stop
            - interpolation = if true, the interpolation by Voight et al. is used, if false, a step function more appropriate for sequence data is used
            -maxDist: the maximal distance to integrate to
        """
        id,haplotypes=self._get_default(singleSite=True, **kwargs)

        if(len(id)==1): id = id[0]
        ind=haplotypes 

        pos=self._seg_sites[id]
        #check if there are more than 2 haplotypes, otherwise return 0
        if len(ind)<2:
            return 0,-1#pos
        #check if we have the topmost snp
        if id==self.nSegsites-1 and dir>0: 
            return 0,self._seg_sites[id]#pos
        if id==0 and dir<0:
            return 0,self._seg_sites[id]

        iEHH=0
        if dir >0:
            dir =1
        else:
            dir = -1

        #if no maxDist is set, set it to +/- infinity
        if maxDist == None:
            maxDist=np.inf
        maxDist = maxDist * dir
        endId               =   self.get_closest_SNP_from_pos(pos+maxDist,direction=-dir)
        curId,curPos,curEHH =   id,self._seg_sites[id],1
        newId,newPos           =   id+dir,self._seg_sites[id+dir]


        while (newId>0 and dir<0) or \
              (newId <self.nSegsites-1 and dir>0):
            delta = (newPos - curPos)*dir
            newEHH = self.EHH(abs(pos-newPos),id=id,onesided=dir,haplotypes=ind)

            if interpolation:
                k = (newEHH+curEHH-.1)*.5*delta
            else:
                k = curEHH*delta #works because curEHH is always >= newEHH

            if verbose:
                print "[pos[%f-%f];id[%i/%i]step: %f/ ehh[i]%f/added: %f/total: %f]"  %(newPos,curPos,newId,curId,delta,newEHH,k,iEHH)
            if newEHH < threshold:
                if interpolation:
                    iEHH        +=  (newPos-curPos) * (1 - (0.05-newEHH)/                       \
                                    (curEHH-newEHH))*.5*(curEHH-0.05)
                    if verbose: print (newPos-curPos)*(1-(0.05-newEHH)/(curEHH-newEHH))*.5*(curEHH-0.05)
                    print iEHH
                    print "end"
                else:
                    iEHH+=curEHH*delta
                break
            iEHH+=k
            #print "[%f][%f-%f][%f/%f]" %(pos,curPos,newPos,curEHH,newEHH)
            curPos,curId,curEHH         =       newPos,newId,newEHH
            newPos,newId                =       self._seg_sites[newId+dir],newId+dir


            #do some adjustment for distance between last snp and end of ehh
            if curId == endId:
                if curId == 0 or curId == self.nSegsites -1:
                    break
                else:
                    #            xmax    -   
                    dist        =   dir * maxDist + dir*(pos-curPos)
                    iEHH += dist*newEHH
                    newPos=(pos+maxDist)
                    break

        print iEHH 
 
        return iEHH,newPos

    def IHS_no_anc_der(self,
                       threshold=.05,verbose=False,
                       interpolation=False,maxDist=None,
                      **kwargs):
        """gets the integrated IHS for all the haplotypes, ignoring
        ancestral/derived alleles"""
        id,haplotypes = self._get_default(singleSite=True,
                                                  **kwargs)
        ihs,dmin                =       self._integrate_EHH(id,haplotypes,-1,threshold,
                                                           maxDist,interpolation,verbose)
        k,dmax                  =       self._integrate_EHH(id,haplotypes,1,threshold,
                                                           maxDist,interpolation,verbose)
        ihs                     +=      k
        return ihs,(dmin,dmax)

    def IHS(self, threshold=0.05,maxDist=None, verbose=False,interpolation=False,noAncDer=False, **kwargs):
        """calculates unstandardized IHS statistic, see Voight et al 2006
            - id = id of the SNP to start at
            - dir = direction, if negative, upstream, if positive, downstream
            - ind = the haplotypes to consider
            - threshold = threshold where to stop
            - interpolation = if true, the interpolation by Voight et al. is used, if false, a step function more appropriate for sequence data is used
            - if it should just be run once
        """

        if noAncDer:
            return self.IHS_no_anc_der(threshold,verbose,interpolation,**kwargs)
        id,haplotypes = self._get_default(singleSite=True, **kwargs)

        iAnc = self.get_haplotypes_with_ancestral_allel(id=id)
        iDer = self.get_haplotypes_with_derived_allel(id=id)

        ihs_derived,dmin = self._integrate_EHH(-1,threshold,maxDist,
                interpolation, verbose,id=id, haplotypes=iDer)
        k,dmax= self._integrate_EHH(1,threshold,maxDist,interpolation,
                                    verbose,id=id, haplotypes=iDer)
        ihs_derived+=k

        ihs_ancestral,amin= self._integrate_EHH(-1,threshold,maxDist,interpolation,verbose, id=id, haplotypes=iAnc)
        k,amax= self._integrate_EHH(1,threshold,maxDist,interpolation,verbose, id=id, haplotypes=iAnc)
        ihs_ancestral+=k

        try:
            ihs = np.log(ihs_ancestral/ihs_derived)
        except ZeroDivisionError:
            ihs = np.inf

        return ihs,ihs_ancestral,ihs_derived,(amin,amax),(dmin,dmax)

    def IHS_external(self):
        np.savetxt("segSites.txt",self._seg_sites,fmt="%f")
        np.savetxt("dump.txt",self._data,fmt="%i")

        s="./ihs3 segSites.txt dump.txt %i > temp.txt" %self.selId
        #print s
        os.system(s)
        try:
            k=np.loadtxt("temp.txt")
        except IOError:
            k=(np.nan,np.nan,np.nan)
        if k[0] ==-1e13:
            k=(np.nan,np.nan,np.nan)
        return (k[2],k[1],k[0])

    def XPEHH(self,x,interpolation=False,verbose=False, **kwargs):
        """untested"""
        id,haplotypes = self._get_default(singleSite=True, **kwargs)

        iPop0 = np.where(self.pops==0)[0]
        iPop1 = np.where(self.pops==1)[0]

        iEHH0,_       =   self._integrate_EHH(id,iPop0,-1,0,x,interpolation,verbose)
        k,_           =   self._integrate_EHH(id,iPop0,1,0,x,interpolation,verbose)
        iEHH0+=k


        iEHH1,_       =   self._integrate_EHH(id,iPop1,-1,0,x,interpolation,verbose)
        k,_           =   self._integrate_EHH(id,iPop1,1,0,x,interpolation,verbose)
        iEHH0+=k

        return(np.log(iEHH0/iEHH1),iEHH0,iEHH1)

#-----------------------------Heterozygosity----------------------------
    def get_SFS(self,**kwargs):
        """avoid double computing global sfs"""
        id,haplotypes=self._get_default(**kwargs)
        if len(id) == self.n_snp() and len(haplotypes) == self.n_hap():
            if self.sfs is None:
                self.sfs = self.make_SFS()
            sfs = self.sfs
        else:
            sfs = self.make_SFS(**kwargs)
        return sfs

    def Heterozygosity(self,**kwargs):
        """calculates heterozygosity directly"""
        sfs = self.get_SFS(**kwargs)
        return sfs.Heterozygosity()

#-----------------------------Pi----------------------------
    def n_pw_differences(self,**kwargs):
        """calculates pi from the SFS"""
        sfs = self.get_SFS(**kwargs)
        return sfs.Pi() 
    Pi = n_pw_differences

#----------------------------Singletons/sfs--------------------------
    def make_SFS(self,**kwargs):
        id,haplotypes=self._get_default(**kwargs)


        nHap=len(haplotypes)
        sfs = SFS(nHap)
        freqs = FreqTable(n_snp=len(id), sample_sizes=[nHap])

        for i in np.array(id):
            snp=self.get_SNP_data(id=i,haplotypes=haplotypes)
            p= sum( snp )
            sfs[p] += 1
            freqs[i] = float(p)

        #sfs=sfs[1:nHap]
        
        return sfs,freqs

    def n_singletons(self,folded=True,**kwargs):
        sfs = self.get_SFS(**kwargs)
        return sfs.n_singletons(folded=folded)

    def hide_singletons(self,folded=True):
        self.hide_allels_with_frequency(freq=1,folded=folded)
    def hide_doubletons(self,folded=True):
        self.hide_allels_with_frequency(freq=2,folded=folded)
    def hide_allels_with_frequency(self,freq=1,folded=True):
        sfs = self.get_SFS(**kwargs)
        toKeep=np.where(np.logical_and(self.freqs!=freq ,
                        np.logical_or(folded, self.nHap-self.freqs!=freq)))[0]
        self._data = self._data[:,toKeep]
        self._seg_sites = np.array(self._seg_sites)[toKeep]
        self.selId = sum(toKeep<self.selId)

        if hasattr(self,"_snp_dict"):
            del self._snp_dict 
        if hasattr(self,"_rec_dict"):
            del self._rec_dict 
            
        self.sfs = None
        self.freqs = None

    def Fay_Wu_H(self,requirePolarizedData=False,**kwargs):
        sfs = self.get_SFS(**kwargs)
        return sfs.Fay_Wu_H(requirePolarizedData=requirePolarizedData)

#-----------------------------S/Tajima's D/Fay+Wu H----------------------------
    def n_segregating_sites(self,**kwargs):
        sfs = self.get_SFS(**kwargs)
        return sfs.n_segregating_sites()
    S = n_segregating_sites

    def Tajimas_D(self,**kwargs):
        """calculates Tajimas D after the formula in Wakeley 2008 p114/115
           checked against the stats program from ms
        """
        sfs = self.get_SFS(**kwargs)
        return sfs.Tajimas_D()
    TD = Tajimas_D


#-----------------------------basic LD stats---------------------------
    def LD_stats(self,ids):
        data=np.empty((self.nHap,2),dtype="bool")

        data[:,0]=self.getDataFromId(ids[0])
        data[:,1]=self.getDataFromId(ids[1])

        n=float(self.nHap)
        p1=sum(data[:,0])/n
        p2=1-p1
        q1=sum(data[:,1])/n
        q2=1-q1
        p11=sum(np.logical_and(data[:,0]==1,data[:,1]==1))/n
        p10=sum(np.logical_and(data[:,0]==1,data[:,1]==0))/n
        p01=sum(np.logical_and(data[:,0]==0,data[:,1]==1))/n
        p00=sum(np.logical_and(data[:,0]==0,data[:,1]==0))/n

        D=p11-p1*q1
        if D<0:
            Dprime=D/np.min(p1*q1,p2*q2)
        elif D>0:
            Dprime=D/np.min(p1*q2,p2*q1)
        else:
            Dprime=0.0

        r=D/np.sqrt(p1*p2*q1*q2)
        return D,Dprime,r
#-----------------------------I/O----------------------------
#----------------------decent quality -----------------------
    def read_next_data_set(self):
        """reads the next data set of a multi dataset mbs output file"""
        if not hasattr(self,"file") or not hasattr(self,"nDataSets"):
            raise ValueError("no DataSet Loaded")

        if self.currDataSet >= self.nDataSets:
            raise ValueError("no more data sets left :-(")

        if not hasattr(self,"type"):
            raise ValueError("I don't know how to handle this input file type")

        if self.type=="mbs":
            self.read_next_mbs_data_set()
        elif self.type=="ms":
            self.read_next_ms_data_set()
        else:
            raise ValueError("Method for this input type not yet implemented")

        self.currDataSet += 1

    def read_next_mbs_data_set(self): 
        """get state of selected Site"
        # note that at selected Site, 1 codes for ancestral, 0 codes for derived
        # allele"""
        h1 = self.file.readline()
        self.selSiteData= [int(x) for x in
                  h1.split(':')[1].replace('d','1').replace('a','0').split()]
        self.nSegsites = int(self.file.readline().split()[1])
        "get pos of segregating sites"
        self._seg_sites = [int(x) for x in self.file.readline().split()[1:]]

        self.selId = sum(np.less(self._seg_sites, self.selPos))
        self._seg_sites.insert(self.selId,self.selPos)
        self._seg_sites = Coords((i,s) for i,s in enumerate(self._seg_sites))

        self._data = np.zeros((self.nHap,self.nSegsites),dtype="bool")
        #self._data=self._data.view(GenArray)
        c=0
        while 1:
            buf=self.file.readline()
            if buf=="\n" or not buf:
                break
            buf = buf.strip()
            s = [int(buf[i]) for i in range(len(buf))]
            s.insert(self.selId,self.selSiteData[c])
            self._data[c,:]=s
            c=c+1

    def read_mbs_file(self,file):
        """reads output from MBS"""
        self.file = open(file)
        self.type = "mbs"

        #first, parse the command as argument
        h1=self.file.readline()
        h1=h1.split()
        tpos,fpos,spos,rpos =(-1,-1,-1,-1)

        for index, value in enumerate(h1):
            if value=='-t': tpos=index
            if value=='-f': fpos=index
            if value=='-s': spos=index
            if value=='-r': rpos=index

        self.parTheta = float(h1[tpos+1])
        self.parRho = float(h1[rpos+1])
        self.nDataSets = int(float(h1[fpos+1]))
        
        self.currDataSet = 0 
        self.selPos = int(float(h1[spos+2]))
        self.length = int(float(h1[spos+1]))
        self.nHap= int(float(h1[2]))
        #self.haplotypes=np.arange(self.nHap)
        self.trajName=h1[spos+3]
        if self.trajName[0:3]=="chr": 
            self.chromosome=self.trajName[3:]
            #print "inferred chr:", self.chromosome


        ln=self.file.readline()
        if ln!="\n":
            print ln
            #assume it is individual ids
            self.individualIDs=ln.split(",")

        #read first DataSet
        self.read_next_data_set()

    def read_ms_file_combine(self,file):
        self.file   =   open(file)
        self.type   =   "msc"
        self.length =   1
        h1          =   self.file.readline()
        h1          =   h1.split()

        self.nHap   =   int(float(h1[1]))
        self.nDataSets= 1
        self.nSegsites= int(float(h1[2]))
        self._data=np.zeros((self.nHap,self.nSegsites),dtype="b")
        #self._data=self._data.view(GenArray)
        self._seg_sites   = np.arange(self.nSegsites,dtype="f")/self.nSegsites
        self._seg_sites = Coords((i,s) for i,s in enumerate(self._seg_sites))

        i = 0 
        while i<self.nSegsites:
            buf=self.file.readline()
            if buf =="":
                break
            if buf[0:8]=="segsites":
                self.file.readline()
                c=0
                while 1:
                    buf=self.file.readline()
                    if buf=="\n" or not buf:
                        break
                    buf = buf.strip()
                    s = int(buf)
                    self._data[c,i]=s
                    c+=1
                i+=1

    def read_ms_file(self,file):
        self.file = open(file)
        self.type = "ms"
        self.length=1

        #parse the command line:
        h1=self.file.readline()
        h1=h1.split()

        self.nHap = int(float(h1[1]))
        self.nDataSets = int(float(h1[2]))
        self.currDataSet = 0
        #self.haplotypes=np.arange(self.nHap)

        self.read_next_data_set()


    def read_next_ms_data_set(self):
        self.reset()
        while 1:
            buf=self.file.readline()
            if buf[0:8]=="segsites":
                buf=buf.split()
                self.nSegsites=int(float(buf[1]))
        #        print self.nSegsites
                break
            if buf == '':
                raise ValueError("could not read input file: number of\
                                     segsites not found")

        #read the positions
        self._seg_sites = [float(x) for x in self.file.readline().split()[1:]]
        #check if all segregating sites are unique, otherwise it will not work:
        if len(np.unique(self._seg_sites)) != len(self._seg_sites):
               raise ValueError("some segregating sites have the same index//ms")
        self._seg_sites = Coords(s for s in self._seg_sites)
        self.coords['pos'] = self._seg_sites
        self.coords['bp']  = self._seg_sites
        self.coords['rec'] = self._seg_sites
        self.coords['id'] = IDCoords()

        #prepare Data
        self._data = np.zeros((self.nHap,self.nSegsites),dtype="int")
        #self._data=self._data.view(GenArray)

        c=0
        while 1:
            buf=self.file.readline()
            if buf=="\n" or not buf:
                break

            buf = buf.strip()
            s = [int(buf[i]) for i in range(len(buf))]
        #    print len(s)
            self._data[c,:]=s
            c=c+1

        #self.set_selected_site(pos=.5)

    def write_to_ms(self,file=None,data=None,segSites=None,par=None,selSiteInData=True):
        """if selSiteInData, the selected Site is writen into the data portion
        as well"""
        if data == None:
            data = self._data
        nHap,nPos = data.shape
        if par == None:
            dummy="dummy"
            if hasattr(self,"chr"): dummy="chr"+self.chr
            par=(nHap,self.parTheta,self.parRho,self.nDataSets,self.length,self.selPos,dummy)
        if segSites == None: 
            segSites = self._seg_sites.values()

        nSNP=len(segSites)


        selSite="".join((str(i)+" " for i in self.selSiteData))
        selSite=selSite.replace("1","d").replace("0","a")

        if file == None:
            #first line fake ms command
            print "./ms %i 1 -t %f -r %f -f 1 %i -s %i %i %s" % par
            print "110 110 110"
            print "\n"
            print "//"
            print "segsites: %i" % nPos
            s="positions: ";
            if segSites != None:
                for ss in segSites:
                    if ss!=self.selPos or selSiteInData:
                        s+=str(ss)+" "
                
            print s 
            for i in range(nHap):
                if selSiteInData:
                    s = "".join([str(int(j)) for j in data[i,:]])
                else:
                    s = "".join([str(int(j)) for j in data[i,:self.selId]])
                    s += "".join([str(int(j)) for j in data[i,self.selId+1:]])
                print s
        else:
            out=open(file,'w')
            #first line fake ms command
            out.write("./ms %i 1 -t %f -r %f -f 1 %i -s %i %i %s\n" %
                      par)
            out.write("110 110 110\n")
            out.write("\n")
            out.write("//\n")
            out.write("segsites: %i\n" % nPos)
            s="positions: ";
            if segSites != None:
                for ss in segSites:
                    if ss!=self.selPos or selSiteInData:
                        s+=str(ss)+" "
                
            out.write(s+"\n") 
            for i in range(nHap):
                if selSiteInData:
                    s = "".join([str(int(j)) for j in data[i,:]])
                else:
                    s = "".join([str(int(j)) for j in data[i,:self.selId]])
                    s += "".join([str(int(j)) for j in data[i,self.selId+1:]])
                out.write(s+"\n") 

    def write_to_mbs(self,file=None,data=None,segSites=None,par=None,selSiteInData=False):
        """if selSiteInData, the selected Site is writen into the data portion
        as well"""
        if data == None:
            data = self._data
        nHap,nPos = data.shape
        if par == None:
            dummy="dummy"
            if hasattr(self,"chr"): dummy="chr"+self.chr
            par=(nHap,self.parTheta,self.parRho,self.nDataSets,self.length,self.selPos,dummy)
        if segSites == None: 
            segSites = self._seg_sites.values()

        nSNP=len(segSites)


        selSite="".join((str(i)+" " for i in self.selSiteData))
        selSite=selSite.replace("1","d").replace("0","a")

        if file == None:
            #first line fake ms command
            print "command:\t./mbs %i -t %f -r %f -f 1 %i -s %i %i %s" % par
            if hasattr(self,"individualIDs"):
                print(",".join(self.individualIDs))
            else:
                print("")
            print "//0-1 allels: "+selSite
            print "segsites: %i" % nPos
            s="positions: ";
            if segSites != None:
                for ss in segSites:
                    if ss!=self.selPos or selSiteInData:
                        s+=str(ss)+" "
                
            print s 
            for i in range(nHap):
                if selSiteInData:
                    s = "".join([str(int(j)) for j in data[i,:]])
                else:
                    s = "".join([str(int(j)) for j in data[i,:self.selId]])
                    s += "".join([str(int(j)) for j in data[i,self.selId+1:]])
                print s
        else:
            out=open(file,'w')
            #first line fake ms command
            out.write("command:\t./mbs %i -t %f -r %f -f 1 %i -s %i %i %s\n" %
                      par)
            if hasattr(self,"individualIDs"):
                out.write(",".join(self.individualIDs))
            else:
                out.write("\n")
            out.write("//0-1 allels: "+selSite+"\n")
            out.write("segsites: %i\n" % nPos)
            s="positions: ";
            if segSites != None:
                for ss in segSites:
                    if ss!=self.selPos or selSiteInData:
                        s+=str(ss)+" "
                
            out.write(s+"\n") 
            for i in range(nHap):
                if selSiteInData:
                    s = "".join([str(int(j)) for j in data[i,:]])
                else:
                    s = "".join([str(int(j)) for j in data[i,:self.selId]])
                    s += "".join([str(int(j)) for j in data[i,self.selId+1:]])
                out.write(s+"\n") 



#--------------------- need improvement ---------------------

    def read_VCF(self, file, panel, selPos, seqPos, 
                    polarizeFile="/data/selectiveSweep/data/ancestral/genes_HC.txt",
                    vcfHasAA=False,excludeNonPolarizable=False):
        """read vcf data into data structure; arguments:
            file:       vcf file name
            panel:      population to use, if it is a list, or np.array
                        they will be used, otherwise determined from header
            selPos:     position of the presumably selected site
            seqPos:     start and ending position
            polarizeFile: file with info on anc/der allele
            vcfHasAA:   if Ancestral Allele is in VCF, this is used instead
            """
        polarize=False
        if not vcfHasAA:
            if polarizeFile != None:
                pf=np.loadtxt(polarizeFile,dtype="S")
                polarize=True
        if type(panel) == list or type(panel) == np.ndarray:
            individualIDs = panel
        else:
            a=np.loadtxt("interim_phase1.20101123.ALL.panel",dtype="S")
            if panel == "ALL":
                individualIDs=[a[i,0]  for i in range(len(a[:,1]))]
            else:
                toKeep=panel
                individualIDs=[a[i,0]  for i in range(len(a[:,1])) if a[i,1] in toKeep]
        print individualIDs

        file= open(file,"r")
        line=file.readline()
        while line[0]=="#":
            header= line
            line=file.readline()

        hd = header.split()
        data = list()
        haplotypesToKeep=[]
        for i,id in enumerate(individualIDs):
            for j,id2 in enumerate(hd):
                if id == id2:
                    haplotypesToKeep.append(j)
            

        #haplotypesToKeep=[i for i in range(len(hd)) if np.array(hd)[i] in individualIDs]  
        print "BLA"
        print haplotypesToKeep
        ls=line.split()
        ht=np.array([ ht[0:3] for ht in np.array(ls)[haplotypesToKeep]])
        snpdata=(np.array([i.split("|") for i in ht])).flatten()
        snpmeta=ls[:9]
        snppos=ls[1]
        snpqual=ls[6]
        snp=(snppos,snpqual,snpdata,snpmeta[3],snpmeta[4])

        data.append(snp)
        nSNPRemoved =0
        while True:
            line=file.readline()
            if not line:
                break
            ls=line.split()
            ht=np.array([ ht[0:3] for ht in np.array(ls)[haplotypesToKeep]])
            snpdata=(np.array([i.split("|") for i in ht])).flatten()
            snpmeta=ls[:9]
            snppos=ls[1]
            snpqual=ls[6]

            if vcfHasAA:
                moreData = ls[7].split( ";" )
                mdd = dict( [l.split("=") for l in moreData] )
                if   mdd["AA"]         == snpmeta[3]: #SNP is ancestral
                    pass
                elif mdd["AA"]         == snpmeta[4]: #SNP is derived
                    snpdata = np.array(1 - np.array(snpdata,dtype="b"),
                                       dtype="S1")
                elif mdd["AA"].upper() == snpmeta[3] and \
                        (not excludeNonPolarizable):
                    pass
                elif mdd["AA"].upper() == snpmeta[4] and \
                        (not excludeNonPolarizable):
                    snpdata = np.array(1 - np.array(snpdata,dtype="b"),
                                       dtype="S1")
                else:
                    if excludeNonPolarizable:
                        nSNPRemoved+=1
                        continue
                    #print snppos, mdd["AA"], snpmeta[3], snpmeta[4]
                    #raise ValueError("SNP neither anc nor derived")
            snp=(snppos,snpqual,snpdata,snpmeta[0],snpmeta[3],snpmeta[4])
            data.append(snp)
        print "removed %d SNP"%nSNPRemoved
        if polarize:
            chr=snpmeta[0][0]
            #refAllele = [i[3] for i in snpmeta] 
            #altAllele = [i[4] for i in snpmeta] 
            #refDict is a dict [chr, pos] => [refH, refC]
            refDict = dict()
            for i,d in enumerate(pf):
                refDict[(pf[i,0],pf[i,1])] = (pf[i,2],pf[i,3])

            toSwitch=[]
            for i,snp in enumerate(data):
                if (snp[3],snp[0]) in refDict:
                    ref = refDict[(snp[3],snp[0])]
                    #print "[%s:%s]1000k: [%s:%s] emilia: [%s:%s]" %(snp[3],snp[0],snp[4],snp[5],ref[0],ref[1])
                    if ref[0] != snp[4]:
                        pass
                        print "snp [%s:%s] is dif between 1k gen and emilia" %  (snp[3],snp[0])
                    elif ref[0] == ref[1]: #snp is polarized correctly
                        if ref[0] == snp[4]:
                            pass
                            print "snp [%s:%s] is polarized correctly" %  (snp[3],snp[0])
                        elif ref[0] == snp[4]:
                            toSwitch.append(snp[0])
                            newData = np.array(["%i"%(1-int(iii)) for iii in data[i][2]],dtype="S1")
                            t = data[i]
                            data[i] = (t[0], t[1],newData,t[3],t[4],t[5])                          
                            print "snp [%s:%s] is polarized wrong!" %  (snp[3],snp[0])
                        else: 
                            pass
                            #print "snp [%s:%s] is ODD (1)!" %  (snp[3],snp[0])
                    else: #snp is not polarized correctly
                        toSwitch.append(snp[0])
                        #print data[i]
                        print "snp [%s:%s] is polarized wrong! (2)" %  (snp[3],snp[0])
                        newData = np.array(["%i"%(1-int(iii)) for iii in data[i][2]],dtype="S1")
                        t = data[i]
                        data[i] = (t[0], t[1],newData,t[3],t[4],t[5])
                        #print data[i]
                else:
                    print "cant polarize %s:%s" % (snp[3],snp[0])

        d2=[d for d in data if (np.any(d[2]=="1") and np.any(d[2]=="0")) or int(d[0])==selPos]
        d3=d2#[d for d in d2 if d[0] == selPos or not np.all(d[2]=="1")]

            
        file.close()    


        nHap=len(d3[0][2])
        nSNP=len(d3)

        segSites=[d[0] for d in d3]
        haplotypes=list()
        for h in range(nHap): haplotypes.append([ d[2][h] for d in d3 ])
        haplotypes=np.array(haplotypes)


        self.type="vcf"
        self.nHap=nHap
        self._seg_sites=np.array([int(i) for i in segSites])
        self.nSegsites=nSNP
        self._data=haplotypes
        #self._data=self._data.view(GenArray)
        #self.haplotypes=np.arange(self.nHap)


        #selId=self.get_closest_SNP_from_pos(selPos)
        try:
            selId=self.getIDfromPos(selPos)
        except:
            raise KeyError(haplotypesToKeep)
        selSiteData=d3[selId][2]
        self.selId=selId
        self.selPos=selPos
        self.selSiteData=selSiteData


        self.parTheta=0
        self.parRho=0
        self.nDataSets=1
        self.length=seqPos[1]-seqPos[0]
        self.selPos-=seqPos[0]
        self._seg_sites=np.array([int(i)-seqPos[0] for i in self._seg_sites])

    def read_IHS(self,file='ihs/sample.ihshap',snp='ihs/sample.ihsmap'):
        """reads input file as used in the ihs program by voight et al., 2006"""
        self.type="ihs"

        data=np.loadtxt(file,dtype="bool")
        snp=np.loadtxt(snp,dtype="S10")

        nHap,nSNP = data.shape
        self.snp=np.array([int(i) for i in snp[:,1]])
        snp=np.array([float(i) for i in snp[:,2]])

        
        self._data=data
        #self._data=self._data.view(GenArray)
        self._seg_sites=snp
        self.nHap=nHap
        self.nSegsites=nSNP

    def read_sweep_phase(self,file='sweep-1.1/CCR5_ceu.phase',
                      snp="sweep-1.1/CCR5_ceu.snp"):
        """reads input file for sweep program in phase format"""

        self.type="phase"
        file=fileinput.input(file)


        data = np.loadtxt(file, dtype="S20")
        data = data[:,2:]
        nHap,nSNP = data.shape
        
        #convert to int
        data = [int(d) for d in data.flatten()]
        data = np.array(data)
        data=data.reshape(nHap,nSNP)

    def read_beagle_phased(self,file='phased/AllDiAllelicSitesEpas1.tibet.like.phased',selPos=0,filter=False):
        """reads a file from the (phased) epas 1 gene, by beagle program"""
        self.type="beagle"

        #get chromosome information
        f=open(file,"r")
        s=f.readline()
        s=f.readline()
        ss=s.split()
        self.chromosome=ss[1].split("_")[0][3:]
        f.close()
        #print "chromosome inferred: ", self.chromosome
        file=fileinput.input(file)

        header = file.readline()
        header = header.split()

        def id2pos(id):
            return int(id.split('_')[1])

        #gets the haplotypes
        for i in range(len(header[3::2])):
            header[2+i*2] += '_a'
            header[3+i*2] += '_b'

        header[1]='pos'

        self.header = header
        self.individualIDs = self.header[2:]

        #annotate samples to pops
        self.han=np.where([s[:2]=="NA" for s in self.individualIDs])[0]
        self.tib =np.where([s[:2]=="DR" or s[:2]=="NQ" for s in
                                  self.individualIDs])[0]
        self.pops=np.zeros(len(self.individualIDs))
        self.pops[self.tib]=1


        type=['S1','I16']; type.extend(['S1' for x in range(len(header)) if x > 1])

        dt={'names':header,'formats':type}
        data = np.loadtxt(file, dtype="S10" ,
                          comments='I', converters = {1:id2pos} )




        #shift to start from 0/x offset
        pos = data[:,1]
        pos = np.asarray([ int(s) for s in pos ])
        #self.minPos = 46358067
        self.minPos = 0
        pos = pos - self.minPos
        self.length = max(pos)

        data= data[:,2:]
        nSNP,nHap = data.shape
        #self.haplotypes=np.arange(nHap)
        data_numeric = np.zeros((nSNP,nHap),dtype='uint8')



        badData = np.zeros(nSNP,dtype="bool")
        self.major = np.empty(nSNP,dtype="S1")
        self.minor = np.empty(nSNP,dtype="S1")
        for i in range(nSNP):
            data_numeric[i,:],badData[i],(self.major[i],self.minor[i]) = makeNumeric(data[i,:])
            if self.verbose: print(i,badData[i])

        file.close()

        #remove all monomorphic SNP
        if filter:
            data_numeric = data_numeric[badData]
            pos = pos[badData]

        self._data = data_numeric.transpose()
        self.nHap,self.nSegsites = self._data.shape

        self._seg_sites = pos
        self.init_SNP_dict()
         #use arbitrary SNP pos for now
        if selPos==0:
            self.setSelectedSite(pos=min(self._seg_sites))
        else:
            self.setSelectedSite(pos=selPos)

    def write_to_haploview(self,file,data=None):
        if data == None:
            data = self._data

        nHap,nPos = data.shape
        families = ['f'+str(i/2) for i in range(nHap)]
        ids = ['i'+str(i/2) for i in range(nHap)]

        output = np.empty((nHap,nPos+2),dtype="S10")
        output[:,0] = np.asarray(families)
        output[:,1] = np.asarray(ids )
        output[:,2:] = data+1


        np.savetxt(file,output,fmt="%s")

    def write_sweep_finder(self,file="sweepfinder_input.txt"):
        output=np.empty((self.nSegsites+1,4),dtype="S10")

    def write_arp(self,file):
        out=open(file,'w')
        out.write( """
        [Profile]
        Title="auto"
        NbSamples=2
        GenotypicData=0
        GameticPhase=1
        DataType=STANDARD
        LocusSeparator=WHITESPACE
        MissingData='?'

        """)
        out.write( """
        [Data]
        [[Samples]]
        SampleName="Tibet"
        """)
        out.write( "SampleSize=%i\n" % (len(self.tib),))
        out.write( "SampleData={\n")

        for i in self.tib:
            out.write( "%i %i %s" %(i, 1, " ".join([str(int(s)) for s in
                                                     self._data[i,:]])))
            out.write("\n")
        out.write( "\n}\n")

        out.write( "SampleName=\"Han\"\n")
        out.write( "SampleSize=%i\n" % (len(self.han),))
        out.write( "SampleData={\n")

        for i in self.han:
            out.write( "%i %i %s" %(i, 1, " ".join([str(int(s)) for s in
                                                     self._data[i,:]])))
            out.write("\n")
        out.write( "\n}\n")

    def write_networkFile(self):
        print "  ;1.0"
        print " ;".join([str(int(s)) for s in self._seg_sites])
        print ";".join(["10" for s in self._seg_sites])

        for i in self.tib:
            print ">tib_%i ;10;0;;;;;;" % i
            print "".join([str(int(s)) for s in self._data[i,:]])
        for i in self.han:
            print ">han_%i ;10;0;;;;;;" % i
            print "".join([str(int(s)) for s in self._data[i,:]])

    def read_arp(self,file="sGeneticsOutput/_spc_samp_1.arp"):
        f       =   open(file)
        line    =   f.readline()
        data    =   []
        sampleId=[]
        while (line !=""):
            if line[0]=="\n" or line[0]=="#":
                line=f.readline()
                continue
            if "SampleName" in line:
                sampleId.append(line.split("\"")[1])
            if "SampleData" in line:
                data.append(f.readline().split()[2:])
                data.append(f.readline().split())
            line=f.readline()
        f.close()

        self._data=np.array(data,dtype="i")
        #self._data=self._data.view(GenArray)
        self.nHap,self.nSegsites=self._data.shape

    def write_sweep(self,file='tibet',chr=2):
        """file for Sabetis sweep program"""

        #write .snp file first, which has 3cols: snpid, chr, pos
        output = np.empty((self.nSegsites+1,3),dtype="S20")
        snpid=["snp_"+str(i) for i in range(self.nSegsites)]
        output[0:] = np.asarray(("snpid","chr","HG16"))
        output[1:,0] = np.asarray(snpid)
        output[1:,1] = np.repeat(chr,self.nSegsites)
        output[1:,2] = np.asarray(self._seg_sites,dtype="S20")
        np.savetxt(file+".snp",output,fmt="%s", delimiter="\t")
        
        #write data file
        output = np.empty((self.nHap,self.nSegsites+2),dtype="S20")
        indid=["ind_"+str(i) for i in range(self.nHap/2)]
        indid=np.repeat(indid,2)
        output[:,0]=indid
        output[:,1]=np.tile(("T","U"),self.nHap/2)
        output[:,2:]=self._data+1
        np.savetxt(file+".phase",output,fmt="%s", delimiter="\t")

        output[0,:]=("position","x","n","folded")
        output[1:,0]=self._seg_sites
        self.make_SFS()
        output[1:,1]=[int(i) for i in self.freqs]
        output[1:,2]=np.repeat(self.nHap,self.nSegsites)
        output[1:,3]=np.repeat(1,self.nSegsites)
        np.savetxt(file,output,fmt="%s",delimiter="\t")

#--------------------- not general ---------------------

    def readNeurospora(self,file="Ellison_2011_SNPdata.txt"):
        self.type="neurospora"
        data=[]
        pos=[]
        chr=[]
        def toInt(c):
            return int(c)
        
        data = np.loadtxt(file,dtype="S20",comments="#", converters = {0:toInt,
                                                                       1:toInt})

        self.chrSNP=data[:,0]
        self._seg_sites=data[:,1]
        data=data[:,2:]
        self.nSegsites,self.nHap=data.shape
        major=np.zeros(self.nSegsites)
        minor=np.zeros(self.nSegsites)
        self._data=np.zeros(data.shape,dtype="i1")
        #self._data=self._data.view(GenArray)
        for i in range(self.nSegsites):
            #self._data[i],b,(major[i],minor[i])    = self.makeNumeric(data[i])
            self._data[i],b,_    = self.makeNumeric(data[i])
            if not b: print i

        self._data=np.transpose(self._data)

    def makeNumeric(self,row):
        nucleotides=np.array(('A','C','G','T'))
        b=1
        alleleFreq=np.array(((sum(row=='A'),sum(row=='C'),sum(row=='G'),sum(row=='T'))))
        if(sum(alleleFreq != 0) != 2):
            if self.verbose: print 'WAAAAAAAAAAAAAAAH: Bad genotype',alleleFreq
            b=0
        k=np.argsort(alleleFreq)[::-1]
        major=nucleotides[k[0]]
        minor=nucleotides[k[1]]
        if alleleFreq[k[1]]==0:
            minor="N"
        nuc=nucleotides[k]
        return (np.where(row==major,0,1),b,(major,minor))



#---------------- ABOVE ARE FUNCTIONS THAT SHOULD BE HERE



        #self.haplotypes=np.arange(self.nHap)













#--- stuff that is questionable ------------------------
    def setSequencedRegionTrack(self,seqFile,format="anna"):
        """this function defines for which region information is available, to
        remove unsequenced region
        file has the following format:
            1. row: header
            2-Nth row: data
            1. col: id (not read)
            2. col type (not read)
            3. col: chromosome (read if chromosome info is available)
            4. col: start of a fragment
            5. col: end of a fragment
        output is saved in self.seqRegion as a tuples of start/endposlist
        """
        """ history: 30.8.2011: created"""
        """ history: 10.10.2011: fixed a bug with readChr"""
        
        #first check if we have chr information
        readChr=False
        if hasattr(self,"chromosome"):
            readChr=True
        
        #read File
        fp = open(seqFile,"r")
        fp.readline() #discard header
        minPos=min(self._seg_sites)
        maxPos=max(self._seg_sites)
        self.seqRegion=([],[])
        for l in fp:
            line = l.split()
            if readChr:
                if self.chromosome != line[3]: continue
            il4 = int(line[4])
            il5 = int(line[5])
            if (il4 > minPos and il4 < maxPos) or \
               (il5 > minPos and il5 < maxPos):
                self.seqRegion[0].append(int(line[4]))
                self.seqRegion[1].append(int(line[5]))

    def removeMonomorphicSites(self):
        self.make_SFS()
        toKeep = np.where(self.freqs!=0)[0]
        self._seg_sites=np.array(self._seg_sites)[toKeep]
        if hasattr(self,"polarizable"): self.polarizable=self.polarizable[toKeep]
        self._data=self._data[:,toKeep]
        self.nSegsites=len(self._seg_sites)
        self.setSelectedSite(pos=self.selPos)

        del self._snp_dict
        self.init_SNP_dict()

    def removeSNP(self,id):
        toKeep=np.ones(self.nSegsites,dtype="Bool")
        toKeep[id]=0
        toKeepId = np.where(toKeep)[0]
        self._seg_sites=np.array(self._seg_sites)[toKeepId]
        if hasattr(self,"polarizable"): self.polarizable=self.polarizable[toKeepId]
        self._data=self._data[:,toKeepId]
        self.nSegsites=len(self._seg_sites)
        self.setSelectedSite(pos=self.selPos)


        del self._snp_dict
        self.init_SNP_dict()
        
    def removeSNPNotSequenced(self):
        """removes the SNP outside of sequenced regions, e.g. because they are
        from a different source. 
        """

        if not hasattr(self,"seqRegion"):
            return

        #get sites that should remain in data
        toKeep=np.zeros(self.nSegsites,dtype="Bool")
        self._seg_sites=np.array(self._seg_sites)
        for i,sp in enumerate(self.seqRegion[0]):
            ep = self.seqRegion[1][i]
            kp = np.logical_and(self._seg_sites>sp ,self._seg_sites<ep)
            toKeep = np.logical_or(kp,toKeep)

        toKeepId = np.where(toKeep)[0]
        self._seg_sites=np.array(self._seg_sites)[toKeepId]
        if hasattr(self,"polarizable"): self.polarizable=self.polarizable[toKeepId]
        self._data=self._data[:,toKeepId]
        self.nSegsites=len(self._seg_sites)
        self.setSelectedSite(pos=self.selPos)


        del self._snp_dict
        self.init_SNP_dict()
        
        return toKeep



#-----------------------------make subset----------------------------
    def get_n_segregating_sitesubset(self,pos=None,id=None,haplotypes=None,setNewSelSite=False):
        """gets a subset of the data as a new mbsData object. if setNewSelSite
        is true, if the selected Site is not in the subset, the next best site
        is put as selected site"""
        pos,id,haplotypes=self._get_default(pos, id, haplotypes)

        if isinstance(id,tuple):
            id=np.arange(id[0],id[1]+1)

        subset = mbsData()
        #stuff that all data sets should have
        subset.data=self._data[haplotypes,:][:,id]
        subset.haplotypes = np.arange(len(haplotypes))
        subset.length=max(self._seg_sites)
        subset.currDataSet=1
        subset.nDataSets=1
        subset.nHap,subset.nSegsites=subset.data.shape
        subset.segSites=np.array(self._seg_sites)[id]
        subset.type=self.type
        subset.verbose=self.verbose

        #selected Site
        if sum(self.selId==id)==1:
            subset.selId=np.where(self.selId==id)[0]
            subset.selPos=self.selPos
            subset.selSiteData=np.array(self.selSiteData)[haplotypes]
        elif setNewSelSite:
            subset.selId=subset.get_closest_SNP_from_pos(self.selPos)
            subset.selPos=subset.get_pos_from_id(subset.selId)
            subset.selSiteData=subset.getDataFromId(subset.selId)

            
        #stuff specific to beagle file
        subset.major=self.major[id]
        subset.minor=self.minor[id]
        subset.individualIDs=np.array(self.individualIDs)[haplotypes]
                               
        subset.han=np.where([s[:2]=="NA" for s in subset.individualIDs])[0]
        subset.tib = tib2=np.where([s[:2]=="DR" or s[:2]=="NQ" for s in
                                                            subset.individualIDs])[0]

        return(subset)

    def getTib(self):
        return self.get_n_segregating_sitesubset(haplotypes=self.tib,setNewSelSite=True)

    def getHan(self):
        return self.get_n_segregating_sitesubset(haplotypes=self.han,setNewSelSite=True)





    def readRefSeqEmilia(self,file):
        rawData=np.loadtxt(file,dtype="S")
        pos         =   np.array([int(i) for i in rawData[:,1]])
        refRaw      =   rawData[:,2]
        for i,d in enumerate(refRaw):
            if d=="NA":
                refRaw[i]=False
        a=dict([(j,refRaw[i])for i,j in enumerate(pos)])
        for i in a:
            if a[i]=="False":
                a[i]=False
        self.ancestralAllel=a


    def readRefSeq(self,file="EPAS1_BGIsnp_orthNuc"):
        rawData=np.loadtxt(file,dtype="S")
        pos         =   np.array([int(i) for i in rawData[1:,1]])
        refHum      =   rawData[1:,2]
        refChimp    =   rawData[1:,3]
        refMaq      =   rawData[1:,4]
        
        polarizable =   refChimp==refMaq
         
        ancestralAllel=dict()
        for i in np.arange(len(pos)):
            if polarizable[i] and refChimp[i] != '-':
                ancestralAllel[pos[i]]=refChimp[i]
            else:
                ancestralAllel[pos[i]]=False
        self.ancestralAllel=ancestralAllel

    def polarizeData(self,keepNA=True,verbose=True):
        if not hasattr(self,"ancestralAllel"):
            raise ValueError("no Data available. Use readRefSeq() to load a file")

        aa=self.ancestralAllel

        toRemove=[]
        self.polarizable=np.zeros(self.nSegsites,dtype="Bool")
        for i in np.arange(len(self._seg_sites)):
            snp=self._seg_sites[i]
            if not snp in aa:
                if verbose: print "SNP %i(%i) has no polarization info:\
                   [m:%s]"\
                %(i,self._seg_sites[i],self.major[i])
                if not keepNA:
                    toRemove.append(i)
                continue
            if aa[snp] == False:
                if verbose: print "SNP %i(%i) is missing data:\
                   [m:%s],[a:%s]"\
                %(i,self._seg_sites[i],self.major[i],aa[snp])
                if sum(self._data[:,i])>self.nHap/2:
                    self._data[:,i] = 1-self._data[:,i]
                    print sum(self._data[:,i])
                if not keepNA:
                    toRemove.append(i)
                continue
            if self.major[i] == aa[snp]:
                if verbose: print "SNP %i(%i) is polarized as expected:\
                   [m:%s],[a:%s]"\
                %(i,self._seg_sites[i],self.major[i],aa[snp])
                self.polarizable[i] = True
            elif self.minor[i] == aa[snp]:
                if verbose: print "SNP %i(%i) is polarized differently:\
                   [m:%s],[a:%s]"\
                %(i,self._seg_sites[i],self.major[i],aa[snp])
                self.polarizable[i] = True
                
                #adjust data
                self._data[:,i] = 1-self._data[:,i]
            else:
                if verbose: print "SNP %i(%i) is different:\
                   [maj:%s],[min:%s],[m[a:%s]"\
                %(i,self._seg_sites[i],self.major[i],self.minor[i],aa[snp])
                if not keepNA:
                    toRemove.append(i)
        for i in toRemove[::-1]:
            self.removeSNP(i)

    def polarizeFromDBSNP(self,f, offset=0):
        """f should be a list of positions of SNP to switch, nothing else
           offset is the starting point if the SNP are in genomic and not local
           coordinates"""
        positions=np.loadtxt(f)
        for i in positions:
            if not i-offset in self._seg_sites:
                continue
            id=self.getIDfromPos(i-offset)
            self._data[:,id] = 1 - self._data[:,id]
            print "SNP " ,id, "/",i,"/",i-offset, " adjusted"

    def simplifyIndividualIDs(self):
        """function to simplify individualIDs to remove the annoying beagle
        attachments"""
        self.individualIDs=[s.split(".")[0] for s in self.individualIDs]
        self.individualIDs=[s.replace("-",".") for s in self.individualIDs]
        self.individualDict=dict([(s,((i-1),i)) for i,s in enumerate(self.individualIDs)])




    def addRecombinationMap(self,startPos=39586954,
                             file="/data/selectiveSweep/IFNL/DATA/genetic_map_GRCh37_chr19.txt"):
        """
        calculates the position of each SNP on a recombination map. This is done
        very lazily in an O(n*m) algorithm, where n is the number of map entries
        and m is the number of SNP. Should be improved when this becomes time
        critical
        """
        geneticMap=[]
        f = open(file)
        for line in f:
            line = line.split()
            curLine=[line[1],line[4]]
            geneticMap.append(curLine)

        geneticMap=np.array(geneticMap,dtype="f4")
        
        segSites_Rec=[]
        for c,curPos_Physical in enumerate(self._seg_sites):
            print c,"/",self.nSegsites
            i=0
            #lazily ignoring the case where the data extends over the chromosome
            while curPos_Physical+startPos > geneticMap[i][0]:
                i+=1
            frac = (curPos_Physical+startPos - geneticMap[i-1][0]) / \
                   (geneticMap[i][0] - geneticMap[i-1][0])
            curPos_Rec = geneticMap[i-1] + frac * (geneticMap[i] - geneticMap[i-1])
            segSites_Rec.append(curPos_Rec[1])




