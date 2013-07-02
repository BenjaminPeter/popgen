#!/usr/bin/env python

"""
    class to handle subsets of haplotypes for sliding window analyses
"""

class HaplotypeSets():
    def setSets(self,sets,includeGlobal=True):
        self.sets=sets
        self.setsIncludeGlobal=includeGlobal

    def setExpandingWindowSets(self,step=1000,startEndPos=None):
        self.windowSize=step
        if startEndPos == None and hasattr(self,"seqRegion"):
            startEndPos=(self.seqRegion[0][0],self.seqRegion[1][-1])
        elif startEndPos == None:
            startEndPos=(min(self.segSites),max(self.segSites))
        startPos=startEndPos[0]
        endPos=startEndPos[1]
        
        startPoints=np.array([startPos for i in range(startPos,endPos,step)])
        endPoints=np.array([i for i in range(startPos,endPos,step)])
        midPoints=(startPoints+endPoints)/2


        self.sets=[]
        #now get SNP in each window
        for i,sp in enumerate(startPoints):
            ep = endPoints[i]
            self.sets.append(np.where(np.logical_and(self.segSites>sp,self.segSites<ep))[0])

        self.windows=(np.array(startPoints),np.array(endPoints),np.array(midPoints))


        self.setsIncludeGlobal=False
        self.createSFS()
        self.createSFSForSets()


    def setSlidingWindowSets(self,windowSize=10000,offset=5000,startEndPos=None,removeEmpty=False):
        """this function calculates sliding windows and saves them in sets.
        if seqRegion is set, it also calculates for each window how much of it is
        covered by each window
        startEndPos is a tuple(startPos,endPos) which gives the start and end of
        the sliding window. if it set to None, it is inferred
        """
        """history: 30.8.2011: created"""

        self.windowSize=windowSize
        if startEndPos == None and hasattr(self,"seqRegion"):
            startEndPos=(self.seqRegion[0][0],self.seqRegion[1][-1])
        elif startEndPos == None:
            startEndPos=(min(self.segSites),max(self.segSites))
        startPos=startEndPos[0]
        endPos=startEndPos[1]
        
        midPoints=np.arange(startPos+windowSize/2.,endPos-windowSize/2.+0.001,offset)
        startPoints=midPoints-windowSize/2
        endPoints=midPoints+windowSize/2


        self.sets=[]
        #now get SNP in each window
        for i,sp in enumerate(startPoints):
            ep = endPoints[i]
            self.sets.append(np.where(np.logical_and(self.segSites>sp,self.segSites<ep))[0])

        #for each window get the coverage of sequence
        coverage=[]
        if hasattr(self,"seqRegion"):
            nSeqRegion=len(self.seqRegion[0])
            c=0 #current seqRegion
            for i,startWin in enumerate(startPoints):
                endWin = endPoints[i]
                cv=0
                for j,startSeq in enumerate(self.seqRegion[0]):
                    endSeq=self.seqRegion[1][j]


                    #sequencing ends  before window or starts after window
                    if endSeq < startWin or startSeq > endWin:
                        #print "result: no overlap",cv
                        continue

                    #if sequence is fully contained in window, add full seq size
                    elif startWin <= startSeq and endSeq <= endWin:
                        cv+=endSeq-startSeq
			print "[%i - %i] vs [%i -%i]" %\
                        (startWin,endWin,startSeq,endSeq),
                        print "result: fully contained",cv
                    elif startSeq <= startWin and endWin <= endSeq:
                        cv+=endWin -startWin
                    #else, right overlap
                    elif startSeq <= startWin and endSeq <= endWin:
                        cv+=endSeq-startWin
			print "[%i - %i] vs [%i -%i]" %\
                        (startWin,endWin,startSeq,endSeq),
                        print "result: right overlap",cv
                    #else, left overlap
                    elif startWin < startSeq and endWin < endSeq:
                        cv+=endWin-startSeq
                        print "[%i - %i] vs [%i -%i]" %\
                        (startWin,endWin,startSeq,endSeq),
                        print "result: left overlap",cv
                    else:
                        pass#print "WAAAAAAAAAAAAAAAAAH", cv
                coverage.append(cv)
        else: coverage = np.array([windowSize for i in midPoints])

        #finally,remove windows with 0 coverage or 0 snp
        toRemove=[]
        for i,s in enumerate(self.sets):
            cv = coverage[i]
            if cv == 0 or len(s) ==0:
                toRemove.append(i)

        midPoints = midPoints.tolist()
        startPoints = startPoints.tolist()
        endPoints = endPoints.tolist()
        if removeEmpty:
            for i in toRemove[::-1]:
                self.sets.pop(i)
                coverage.pop(i)
                midPoints.pop(i)
                startPoints.pop(i)
                endPoints.pop(i)

        self.coverage=np.array(coverage)
        self.windows=(np.array(startPoints),np.array(endPoints),np.array(midPoints))


        self.setsIncludeGlobal=False
        self.createSFS()
        self.createSFSForSets()
    def getEHHforSetsAvg(self,x):
        """"""
        #first, calculate all EHH values
        ehh=np.empty(self.nSegsites)
        for i,ss in enumerate(self.segSites):
            ehh[i]=self.getEHH(x,id=i,derivedOnly=True)

        #then get avg for sets
        ehhSets=np.empty(len(self.sets))
        for i,set in enumerate(self.sets):
            ehhSets[i]=np.mean(ehh[set])
        self.ehh=ehhSets
        return self.ehh

    def getIHSforSetsAvg(self,requirePolarizedData=False):
        ihs=np.empty(self.nSegsites)
        for i,ss in enumerate(self.segSites):
            ihs[i]=self.getIHS(id=i)[0]
        ihs=np.ma.array(ihs)
        ihs[np.logical_or(ihs==np.inf,-ihs==np.inf)]=np.ma.masked
        if requirePolarizedData:
            ihs[np.logical_not(self.polarizable)] = np.ma.masked


        ihsSets=np.empty(len(self.sets))
        for i,set in enumerate(self.sets):
            
            ihsSets[i]=np.mean(ihs[set])
        self.ihs=ihsSets
        return self.ihs

    def getHFromSets(self):
        if not hasattr(self,"pi"):
            self.getPiFromSets()

        if not hasattr(self,"S"):
            self.getSFromSets()

        return self.pi/self.S
            
    def get_PiFromSets(self):
        if self.setsIncludeGlobal:
            self.pi=np.empty(len(self.sets)+1)
        else:
            self.pi=np.empty(len(self.sets))

        n=self.nHap
        denominator = n*(n-1.0)

        if self.setsIncludeGlobal:
            self.pi=np.empty(len(self.sets)+1)
        else:
            self.pi=np.empty(len(self.sets))

        for k,sfs in enumerate(self.setSFS):
            thetaPi=0.0
            for i,freq in enumerate(sfs):
                if i==0 or i==self.nHap:
                    continue
                thetaPi+= 2.0*i*(n-i)*freq
            self.pi[k] = thetaPi/denominator
        return self.pi

    def createSFSForSets(self):
        self.setSFS=[]
        for s in self.sets:
            condFreq=self.freqs[s]
            sfs=np.zeros(self.nHap+1)
            for i in condFreq:
                sfs[i]+=1
            self.setSFS.append(sfs)
        if self.setsIncludeGlobal:
            self.setSFS.insert(0,self.sfs)
        self.setSFS=np.array(self.setSFS)

    def getNSingletonsFromSets(self):
        return self.setSFS[:,1]

    def getFayWuHFromSets(self,requirePolarizedData=False):
        if requirePolarizedData:
            #copy sets
            tempSets=[s.copy() for s in self.sets]

            #polarizableId=np.where(self.polarizable)[0]
            for i,s in enumerate(self.sets):
                self.sets[i]=np.ma.array(s)
                for j,snp in enumerate(s):
                    if not self.polarizable[snp]:
                        self.sets[i][j]=np.ma.masked
            self.sets=[x.compressed() for x in self.sets]
            if hasattr(self,"pi"):
                del self.pi
            if hasattr(self,"setSFS"):
                del self.setSFS
            if hasattr(self,"freqs"):
                del self.freqs
            if hasattr(self,"sfs"):
                del self.sfs
            self.make_SFS()
            self.make_SFSForSets()
            result = self.getFayWuHFromSets()

            #no undo all the stuff
            del self.pi
            del self.setSFS
            del self.freqs
            del self.sfs
            self.sets=tempSets
            self.make_SFS()
            self.make_SFSForSets()
            return result
        else:
            if not hasattr(self,"pi"):
                self.get_PiFromSets()
            thetaH=0.0
            n=self.nHap
            denominator = n*(n-1.0)

            if self.setsIncludeGlobal:
                self.fwh=np.empty(len(self.sets)+1)
            else:
                self.fwh=np.empty(len(self.sets))

            for k,sfs in enumerate(self.setSFS):
                thetaH=0.0
                for i,freq in enumerate(sfs):
                    if i==0 or i==self.nHap:
                        continue
                    thetaH+= 2*i*i*freq

                self.fwh[k] = self.pi[k] - thetaH/denominator
            return self.fwh

    def getSFromSets(self):
        if self.setsIncludeGlobal:
            self.S=np.empty(len(self.sets)+1,dtype="I")
        else:
            self.S=np.empty(len(self.sets),dtype="I")
        
        for k,sfs in enumerate(self.setSFS):
            self.S[k] = sum(self.setSFS[k][1:]) 
        return self.S

    def getTajimasDfromSets(self):
        """ """
        n=self.nHap

        if not hasattr(self,"pi"):
            self.get_PiFromSets()

        if not hasattr(self,"S"):
            self.get_N_segregating_sitesFromSets()

        if self.setsIncludeGlobal:
            self.td=np.empty(len(self.sets)+1)
        else:
            self.td=np.empty(len(self.sets))

        #coefficient a1
        a1=0.0
        for i in range(1,n):
            a1+=1./i
        #coefficient a2
        a2=0.
        for i in range(1,n):
            a2+=1./(i*i)

        #more coeffs
        e1 = 1./a1 * ((n+1.)/(3.*(n-1.))-1/a1)
        e2 = 1./(a1*a1+a2) * (((2.*(n*n+n-3.))/(9.*n*(n-1))) - ((n+2.)/(n*a1))+a2/(a1*a1))
        #done with coeff
        
        for k,sfs in enumerate(self.setSFS):
            tdRaw=self.pi[k] - self.S[k]/a1
            vartD=e1*self.S[k]+e2*self.S[k]*(self.S[k]-1.0)
            self.td[k]=tdRaw/np.sqrt(vartD)

        return self.td
