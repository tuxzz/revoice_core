import ctypes
import numpy as np
import numpy.ctypeslib as npct
from . import rtpyin

dll = ctypes.CDLL("librevoice.dll")
RvReal = ctypes.c_double
RvReal_1d = npct.ndpointer(dtype = np.float64, ndim = 1, flags = "C")
RvReal_2d = npct.ndpointer(dtype = np.float64, ndim = 2, flags = "C")

class RvRTMonoPitchProcessorParameter(ctypes.Structure):
    _fields_ = [
        ("samprate", RvReal),
        ("maxTransSemitone", RvReal), ("maxminFreq", RvReal),
        ("transSelf", RvReal), ("yinTrust", RvReal), ("energyThreshold", RvReal),
        ("hopSize", ctypes.c_int), ("nSemitone", ctypes.c_int),
        ("binPerSemitone", ctypes.c_int), ("maxObsLength", ctypes.c_int),
    ]

class RvRTMonoPitchProcessor(ctypes.Structure):
    pass

pRvRTMonoPitchProcessorParameter = ctypes.POINTER(RvRTMonoPitchProcessorParameter)
pRvRTMonoPitchProcessor = ctypes.POINTER(RvRTMonoPitchProcessor)

rvCreateRTMonoPitchProcessorParameter = dll.rvCreateRTMonoPitchProcessorParameter
rvCreateRTMonoPitchProcessorParameter.argtypes = [ctypes.c_int, RvReal, ctypes.c_int, RvReal, RvReal]
rvCreateRTMonoPitchProcessorParameter.restype = pRvRTMonoPitchProcessorParameter

rvCreateRTMonoPitchProcessorParameterFromRTPYin = dll.rvCreateRTMonoPitchProcessorParameterFromRTPYin
rvCreateRTMonoPitchProcessorParameterFromRTPYin.argtypes = [rtpyin.pRvRTPYinProcessorParameter]
rvCreateRTMonoPitchProcessorParameterFromRTPYin.restype = pRvRTMonoPitchProcessorParameter

rvDestroyRTMonoPitchProcessorParameter = dll.rvDestroyRTMonoPitchProcessorParameter
rvDestroyRTMonoPitchProcessorParameter.argtypes = [pRvRTMonoPitchProcessorParameter]
rvDestroyRTMonoPitchProcessorParameter.restype = None

rvCreateRTMonoPitchProcessor = dll.rvCreateRTMonoPitchProcessor
rvCreateRTMonoPitchProcessor.argtypes = [pRvRTMonoPitchProcessorParameter]
rvCreateRTMonoPitchProcessor.restype = pRvRTMonoPitchProcessor

rvRTMonoPitchParam = dll.rvRTMonoPitchParam
rvRTMonoPitchParam.argtypes = [pRvRTMonoPitchProcessor]
rvRTMonoPitchParam.restype = pRvRTMonoPitchProcessorParameter

rvCallRTMonoPitch = dll.rvCallRTMonoPitch
rvCallRTMonoPitch.argtypes = [pRvRTMonoPitchProcessor, RvReal_1d, RvReal_2d, ctypes.c_int, RvReal_1d]
rvCallRTMonoPitch.restype = ctypes.c_int

rvMonoPitchNextOutputLength = dll.rvMonoPitchNextOutputLength
rvMonoPitchNextOutputLength.argtypes = [pRvRTMonoPitchProcessor]
rvMonoPitchNextOutputLength.restype = ctypes.c_int

rvDestroyRTMonoPitchProcessor = dll.rvDestroyRTMonoPitchProcessor
rvDestroyRTMonoPitchProcessor.argtypes = [pRvRTMonoPitchProcessor]
rvDestroyRTMonoPitchProcessor.restype = None

def parameterFromPYin(pyin):
    hopSize = pyin.hopSize
    samprate = pyin.samprate
    nSemitone = int(np.ceil(np.log2(pyin.maxFreq / pyin.minFreq) * 12.0))
    maxTransSemitone = (pyin.hopSize / pyin.samprate) / (256.0 / 44100.0) * 3.0
    minFreq = pyin.minFreq
    return hopSize, samprate, nSemitone, maxTransSemitone, minFreq

class Processor:
    def __init__(self, hopSize, samprate, nSemitone, maxTransSemitone, minFreq, **kwargs):
        self.maxObsLength = kwargs.get("maxObsLength", 128)

        self.hopSize = int(hopSize)
        self.samprate = samprate
        self.nSemitone = nSemitone
        self.maxTransSemitone = maxTransSemitone
        self.minFreq = minFreq
        self.binPerSemitone = kwargs.get("binPerSemitone", 5)
        self.transSelf = kwargs.get("transSelf", 0.999)
        self.yinTrust = kwargs.get("yinTrust", 0.5)
        self.energyThreshold = kwargs.get("energyThreshold", 1e-8)
        
        param = rvCreateRTMonoPitchProcessorParameter(hopSize, samprate, nSemitone, maxTransSemitone, minFreq)
        param.contents.binPerSemitone = self.binPerSemitone
        param.contents.transSelf = self.transSelf
        param.contents.yinTrust = self.yinTrust
        param.contents.energyThreshold = self.energyThreshold
        param.contents.maxObsLength = self.maxObsLength
        self.proc = rvCreateRTMonoPitchProcessor(param)
        rvDestroyRTMonoPitchProcessorParameter(param)
    
    def __del__(self):
        rvDestroyRTMonoPitchProcessor(self.proc)

    def __call__(self, x, obsProb):
        if(len(x) != self.hopSize * 2):
            raise ValueError("length of x must be hopSize * 2")
        if(obsProb.ndim != 2 or obsProb.shape[1] != 2):
            raise ValueError("invalid obsProb")
        n = rvMonoPitchNextOutputLength(self.proc)
        o = np.zeros(n, dtype = np.float64)
        realN = rvCallRTMonoPitch(self.proc, x, obsProb, obsProb.shape[0], o)
        '''
        if(obsProb.shape[0] == 0):
            v = rvCallRTMonoPitch.argtypes
            rvCallRTMonoPitch.argtypes = [pRvRTMonoPitchProcessor, RvReal_1d, ctypes.POINTER(RvReal), ctypes.c_int, RvReal_1d]
            realN = rvCallRTMonoPitch(self.proc, x, None, 0, o)
            rvCallRTMonoPitch.argtypes = v
        else:
            realN = rvCallRTMonoPitch(self.proc, x, obsProb, obsProb.shape[0], o)'''
        assert n == realN
        
        return o