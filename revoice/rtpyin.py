import ctypes
import numpy as np
import numpy.ctypeslib as npct

from .common import *

dll = ctypes.CDLL("librevoice.dll")
RvReal = ctypes.c_double
RvReal_1d = npct.ndpointer(dtype = np.float64, ndim = 1, flags = "C")
RvReal_2d = npct.ndpointer(dtype = np.float64, ndim = 2, flags = "C")

class RvRTPYinProcessorParameter(ctypes.Structure):
    _fields_ = [
        ("pdf", RvReal_1d),
        ("samprate", RvReal),
        ("minFreq", RvReal), ("maxFreq", RvReal),
        ("valleyThreshold", RvReal), ("valleyStep", RvReal),
        ("probThreshold", RvReal), ("weightPrior", RvReal),
        ("bias", RvReal),
        ("hopSize", ctypes.c_int), ("maxWindowSize", ctypes.c_int),
        ("maxIter", ctypes.c_int), ("pdfSize", ctypes.c_int),
        ("maxInputSegment", ctypes.c_int),
        ("prefilter", ctypes.c_bool), ("isPdfDefault", ctypes.c_bool),
    ]

class RvRTPYinProcessor(ctypes.Structure):
    pass

pRvRTPYinProcessorParameter = ctypes.POINTER(RvRTPYinProcessorParameter)
pRvRTPYinProcessor = ctypes.POINTER(RvRTPYinProcessor)

rvCreateRTPYinProcessorParameter = dll.rvCreateRTPYinProcessorParameter
rvCreateRTPYinProcessorParameter.argtypes = [RvReal, RvReal, RvReal, RvReal_1d, ctypes.c_int]
rvCreateRTPYinProcessorParameter.restype = pRvRTPYinProcessorParameter

rvDestroyRTPYinProcessorParameter = dll.rvDestroyRTPYinProcessorParameter
rvDestroyRTPYinProcessorParameter.argtypes = [pRvRTPYinProcessorParameter]
rvDestroyRTPYinProcessorParameter.restype = None

rvRTPYinBufferUsed = dll.rvRTPYinBufferUsed
rvRTPYinBufferUsed.argtypes = [pRvRTPYinProcessor]
rvRTPYinBufferUsed.restype = ctypes.c_int

rvRTPYinDumpBuffer = dll.rvRTPYinDumpBuffer
rvRTPYinDumpBuffer.argtypes = [pRvRTPYinProcessor, RvReal_1d]
rvRTPYinDumpBuffer.restype = None

rvCreateRTPYinProcessor = dll.rvCreateRTPYinProcessor
rvCreateRTPYinProcessor.argtypes = [pRvRTPYinProcessorParameter]
rvCreateRTPYinProcessor.restype = pRvRTPYinProcessor

rvRTPYinParam = dll.rvRTPYinParam
rvRTPYinParam.argtypes = [pRvRTPYinProcessor]
rvRTPYinParam.restype = pRvRTPYinProcessorParameter

rvRTPYinDelayed = dll.rvRTPYinDelayed
rvRTPYinDelayed.argtypes = [pRvRTPYinProcessor]
rvRTPYinDelayed.restype = ctypes.c_int

rvCallRTPYin = dll.rvCallRTPYin
rvCallRTPYin.argtypes = [pRvRTPYinProcessor, RvReal_1d, ctypes.c_int, RvReal_2d, ctypes.c_int]
rvCallRTPYin.restype = ctypes.c_int

rvDestroyRTPYinProcessor = dll.rvDestroyRTPYinProcessor
rvDestroyRTPYinProcessor.argtypes = [pRvRTPYinProcessor]
rvDestroyRTPYinProcessor.restype = None

rvRTPYinDelay = dll.rvRTPYinDelay
rvRTPYinDelay.argtypes = [pRvRTPYinProcessorParameter]
rvRTPYinDelay.restype = ctypes.c_int

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

        self.pdf = kwargs.get("pdf", None)

        self.maxInputSegment = kwargs.get("maxInputSegment", self.hopSize)
        self.maxWindowSize = max(roundUpToPowerOf2(self.samprate / self.minFreq * 4), self.hopSize)
        if(self.pdf is None):
            v = rvCreateRTPYinProcessorParameter.argtypes
            rvCreateRTPYinProcessorParameter.argtypes = [RvReal, RvReal, RvReal, ctypes.POINTER(RvReal), ctypes.c_int]
            param = rvCreateRTPYinProcessorParameter(self.minFreq, self.maxFreq, self.samprate, None, 0)
            rvCreateRTPYinProcessorParameter.argtypes = v
        else:
            param = rvCreateRTPYinProcessorParameter(self.minFreq, self.maxFreq, self.samprate, self.pdf, self.pdf.shape[0])
        param.contents.hopSize = self.hopSize
        param.contents.maxIter = self.maxIter
        param.contents.prefilter = self.prefilter
        param.contents.valleyThreshold = self.valleyThreshold
        param.contents.valleyStep = self.valleyStep
        param.contents.probThreshold = self.probThreshold
        param.contents.weightPrior = self.weightPrior
        param.contents.bias = self.bias
        param.contents.maxInputSegment = self.maxInputSegment
        param.contents.maxWindowSize = self.maxWindowSize
        self.proc = rvCreateRTPYinProcessor(param)
        rvDestroyRTPYinProcessorParameter(param)

        if(self.proc is None):
            raise RuntimeError("cannot create processor")
    
    def __del__(self):
        rvDestroyRTPYinProcessor(self.proc)

    @property
    def delayed(self):
        return rvRTPYinDelay(rvRTPYinParam(self.proc))

    @property
    def buffer(self):
        n = rvRTPYinBufferUsed(self.proc)
        out = np.zeros(n, dtype = np.float64)
        rvRTPYinDumpBuffer(self.proc, out)
        return out

    def __call__(self, x):
        maxOut = 128
        freqProb = np.zeros((maxOut, 2), dtype = np.float64)
        if(x is None):
            v = rvCallRTPYin.argtypes
            rvCallRTPYin.argtypes = [pRvRTPYinProcessor, ctypes.POINTER(RvReal), ctypes.c_int, RvReal_2d, ctypes.c_int]
            nOut = rvCallRTPYin(self.proc, None, 0, freqProb, maxOut)
            rvCallRTPYin.argtypes = v
        else:
            nOut = rvCallRTPYin(self.proc, x, x.shape[0], freqProb, maxOut)
        
        if(nOut == -1):
            return None
        return freqProb[:nOut].copy()