import numpy as np
import pylab as pl

from .common import *
from . import yin, pyin, rtfilter

class Processor:
    def __init__(self, sr, **kwargs):
        self.samprate = float(sr)
        self.hopSize = kwargs.get("hopSize", roundUpToPowerOf2(self.samprate * 0.0025))

        self.minFreq = kwargs.get("minFreq", 80.0)
        self.maxFreq = kwargs.get("maxFreq", 1000.0)
        self.maxIter = kwargs.get("maxIter", 4)
        self.prefilter = kwargs.get("prefilter", True)

        self.valleyThreshold = kwargs.get("valleyThreshold", 1.0)
        self.valleyStep = kwargs.get("valleyStep", 0.01)

        self.probThreshold = kwargs.get("probThreshold", 0.02)
        self.weightPrior = kwargs.get("weightPrior", 5.0)
        self.bias = kwargs.get("bias", 1.0)

        self.pdf = kwargs.get("pdf", pyin.normalized_pdf(1.7, 6.8, 0.0, 1.0, 128))

        self.maxInputSegment = kwargs.get("maxInputSegment", self.hopSize)
        self.maxWindowSize = max(roundUpToPowerOf2(self.samprate / self.minFreq * 4), self.hopSize)
        halfWindowSize = self.maxWindowSize // 2
        self.delay = halfWindowSize - 1
        if(self.prefilter):
            filterOrder = int(2048 * sr / 44100.0)
            if(filterOrder % 2 == 0):
                filterOrder += 1
            kernel = sp.firwin(filterOrder, max(1500.0, self.maxFreq * 4.0), window = "blackman", nyq = sr / 2.0)
            self.prefilterProc = rtfilter.Procressor(kernel)
            self.delay += filterOrder
        self.buffer = np.zeros(halfWindowSize, dtype = np.float64)
        self._internalDelayed = 0
    
    @property
    def delayed(self):
        if(self.prefilter):
            return self._internalDelayed + self.prefilterProc.delayed
        else:
            return self._internalDelayed

    def __call__(self, x):
        halfMaxWindowSize = self.maxWindowSize // 2
        if(x is not None and len(x) > self.maxInputSegment):
            raise ValueError("length of x cannot be greater than maxInputSegment(got %d)" % len(x))

        if(self.prefilter):
            x = self.prefilterProc(x)
        
        if(x is None):
            if(self._internalDelayed > 0):
                self.buffer = np.concatenate((self.buffer, np.zeros(self.hopSize)))
            else:
                return None
        else:
            nX = len(x)
            self._internalDelayed += nX
            self.buffer = np.concatenate((self.buffer, x))
            if(len(self.buffer) < self.maxWindowSize):
                return None
        
        assert len(self.buffer) >= self.maxWindowSize

        windowSize = 0
        newWindowSize = max(roundUpToPowerOf2(self.samprate / self.minFreq * 4), self.hopSize * 2)
        iIter = 0
        while(newWindowSize != windowSize and iIter < self.maxIter):
            windowSize = newWindowSize
            halfDelta = (self.maxWindowSize - windowSize) // 2
            frame = self.buffer[halfDelta:self.maxWindowSize - halfDelta]

            buff = yin.difference(frame)
            buff = yin.cumulativeDifference(buff)
            valleyIndexList = yin.findValleys(buff, self.minFreq, self.maxFreq, self.samprate, threshold = self.valleyThreshold, step = self.valleyStep)
            nValley = len(valleyIndexList)
            if(valleyIndexList):
                possibleFreq = min(self.maxFreq, max(self.samprate / valleyIndexList[-1] - 20.0, self.minFreq))
                newWindowSize = max(int(np.ceil(self.samprate / possibleFreq * 4)), self.hopSize * 2)
                if(newWindowSize % 2 != 0):
                    newWindowSize += 1
                iIter += 1
        
        pdfSize = len(self.pdf)
        freqProb = np.zeros((nValley, 2), dtype = np.float64)
        probTotal = 0.0
        weightedProbTotal = 0.0
        for iValley, valley in enumerate(valleyIndexList):
            ipledIdx, ipledVal = parabolicInterpolation(buff, valley)
            freq = self.samprate / ipledIdx
            v0 = 1 if(iValley == 0) else min(1.0, buff[valleyIndexList[iValley - 1]] + 1e-10)
            v1 = 0 if(iValley == nValley - 1) else max(0.0, buff[valleyIndexList[iValley + 1]]) + 1e-10
            prob = 0.0
            for i in range(int(v1 * pdfSize), int(v0 * pdfSize)):
                prob += self.pdf[i] * (1.0 if(ipledVal < i / pdfSize) else 0.01)
            prob = min(prob, 0.99)
            prob *= self.bias
            probTotal += prob
            if(ipledVal < self.probThreshold):
                prob *= self.weightPrior
            weightedProbTotal += prob
            freqProb[iValley] = freq, prob

        # renormalize
        if(nValley > 0 and weightedProbTotal != 0.0):
            freqProb.T[1] *= probTotal / weightedProbTotal

        self.buffer = self.buffer[self.hopSize:].copy()
        self._internalDelayed -= self.hopSize

        return freqProb