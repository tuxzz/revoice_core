import numpy as np

from . import yin
from . import rtfilter
from .common import *

class Processor:
    def __init__(self, sr, **kwargs):
        self.samprate = float(sr)
        self.hopSize = kwargs.get("hopSize", roundUpToPowerOf2(self.samprate * 0.0025))

        self.minFreq = kwargs.get("minFreq", 80.0)
        self.maxFreq = kwargs.get("maxFreq", 1000.0)
        self.windowSize = kwargs.get("windowSize", max(roundUpToPowerOf2(self.samprate / self.minFreq * 2), self.hopSize * 4))
        self.prefilter = kwargs.get("prefilter", True)
        self.maxInputSegment = kwargs.get("maxInputSegment", self.hopSize)

        self.valleyThreshold = kwargs.get("valleyThreshold", 0.5)
        self.valleyStep = kwargs.get("valleyStep", 0.01)

        halfWindowSize = self.windowSize // 2
        self.delay = halfWindowSize - 1
        if(self.prefilter):
            filterOrder = int(2048 * sr / 44100.0)
            if(filterOrder % 2 == 0):
                filterOrder += 1
            kernel = sp.firwin(filterOrder, max(1250.0, self.maxFreq * 1.25), window = "blackman", nyq = sr / 2.0)
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
        halfWindowSize = self.windowSize // 2
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
            if(len(self.buffer) < self.windowSize):
                return None
        
        assert len(self.buffer) >= self.windowSize

        frame = self.buffer[:self.windowSize]
        buff = yin.difference(frame)
        buff = yin.cumulativeDifference(buff)
        valleyIndexList = yin.findValleys(buff, self.minFreq, self.maxFreq, self.samprate, threshold = self.valleyThreshold, step = self.valleyStep)
        out = self.samprate / parabolicInterpolation(buff, valleyIndexList[-1], val = False) if(valleyIndexList) else 0.0

        self.buffer = self.buffer[self.hopSize:].copy()
        self._internalDelayed -= self.hopSize
        return out