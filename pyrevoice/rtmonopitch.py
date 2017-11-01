import numpy as np

from .common import *
from . import hmm, rthmm, monopitch

class Processor:
    def __init__(self, hopSize, samprate, nSemitone, maxTransSemitone, minFreq, **kwargs):
        self.proc = monopitch.Processor(hopSize, samprate, nSemitone, maxTransSemitone, minFreq, **kwargs)
        self.maxObsLength = kwargs.get("maxObsLength", 128)

        self.hopSize = self.proc.hopSize
        self.samprate = self.proc.samprate
        self.nSemitone = self.proc.nSemitone
        self.maxTransSemitone = self.proc.maxTransSemitone
        self.minFreq = self.proc.minFreq
        self.binPerSemitone = self.proc.binPerSemitone
        self.transSelf = self.proc.transSelf
        self.yinTrust = self.proc.yinTrust
        self.energyThreshold = self.proc.energyThreshold
        self.model = rthmm.SparseHMM(self.proc.model.init, self.proc.model.frm, self.proc.model.to, self.proc.model.transProb, self.maxObsLength)

        self.currObsLength = 0
        self.obsProbList = []
        self.obsSilentList = []

    def __call__(self, x, obsProb):
        if(len(x) != self.proc.hopSize * 2):
            raise ValueError("length of x must be 2 * hopSize")
        
        # constant
        nBin = int(self.nSemitone * self.binPerSemitone)
        nState = len(self.model.init)
        maxFreq = self.minFreq * np.power(2, self.nSemitone / 12)
        
        # feed and decode
        obs = self.proc.calcStateProb(obsProb)
        self.model.feed(obs)
        self.currObsLength = min(self.currObsLength + 1, self.maxObsLength)
        path = self.model.viterbiDecode(self.currObsLength)

        # save state
        frame = simpleDCRemove(x)
        meanEnergy = np.mean(frame ** 2)
        isSilent = meanEnergy < self.energyThreshold

        self.obsProbList.append(obsProb)
        self.obsSilentList.append(isSilent)
        while(len(self.obsProbList) > self.maxObsLength):
            del self.obsProbList[0]
            del self.obsSilentList[0]

        # extract frequency from path
        out = np.zeros(self.currObsLength, dtype = np.float64)
        for iHop in range(self.currObsLength):
            if(path[iHop] < nBin):
                hmmFreq = self.minFreq * np.power(2, path[iHop] / (12.0 * self.binPerSemitone))
                if(len(self.obsProbList[iHop]) == 0):
                    bestFreq = hmmFreq
                else:
                    iNearest = np.argmin(np.abs(self.obsProbList[iHop].T[0] - hmmFreq))
                    bestFreq = self.obsProbList[iHop][iNearest][0]
                    if(bestFreq < self.minFreq or bestFreq > maxFreq or abs(np.log2(bestFreq / self.minFreq) * 12 * self.binPerSemitone - path[iHop]) > 1.0):
                        bestFreq = hmmFreq
            else:
                bestFreq = -self.minFreq * np.power(2, (path[iHop] - nBin) / (12 * self.binPerSemitone))
            out[iHop] = bestFreq
        
        # mark unvoiced->voiced bound as voiced
        for iHop in range(1, self.currObsLength):
            if(out[iHop - 1] <= 0.0 and out[iHop] > 0.0):
                windowSize = max(int(np.ceil(self.samprate / out[iHop] * 4)), self.hopSize * 2)
                if(windowSize % 2 == 1):
                    windowSize += 1
                frameOffset = int(round(windowSize / self.hopSize / 2))
                out[max(0, iHop - frameOffset):iHop] = out[iHop]
        
        # mark silent frame as unvoiced
        for iHop, isSilent in enumerate(self.obsSilentList):
            if(out[iHop] > 0.0 and isSilent):
                out[iHop] = 0.0
        return out