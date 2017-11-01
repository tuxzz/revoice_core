import ctypes
import numpy as np
import numpy.ctypeslib as npct

from .common import *

dll = ctypes.CDLL("librevoice.dll")
RvReal = ctypes.c_double
RvReal_1d = npct.ndpointer(dtype = np.float64, ndim = 1, flags = "C")
RvReal_2d = npct.ndpointer(dtype = np.float64, ndim = 2, flags = "C")

rvFirwinSingleBand = dll.rvFirwinSingleBand
rvFirwinSingleBand.argtypes = [ctypes.c_int, RvReal, RvReal, ctypes.c_char_p, ctypes.c_bool, RvReal, RvReal_1d]
rvFirwinSingleBand.restype = None

def firwinSingleBand(numtaps, passBegin, passEnd, window, nyq, scale = True):
    out = np.zeros(numtaps, dtype = np.float64)
    rvFirwinSingleBand(numtaps, passBegin, passEnd, window.encode("utf-8"), scale, nyq, out)
    return out

class RvRTFilter(ctypes.Structure):
    pass

pRvRTFilter = ctypes.POINTER(RvRTFilter)

rvCreateRTFilter = dll.rvCreateRTFilter
rvCreateRTFilter.argtypes = [RvReal_1d, ctypes.c_int, ctypes.c_int]
rvCreateRTFilter.restype = pRvRTFilter

rvRTFilterNextOutputSize = dll.rvRTFilterNextOutputSize
rvRTFilterNextOutputSize.argtypes = [pRvRTFilter, ctypes.c_int]
rvRTFilterNextOutputSize.restype = ctypes.c_int

rvCallRTFilter = dll.rvCallRTFilter
rvCallRTFilter.argtypes = [pRvRTFilter, RvReal_1d, ctypes.c_int, RvReal_1d]
rvCallRTFilter.restype = ctypes.c_int

rvDestroyRTFilter = dll.rvDestroyRTFilter
rvDestroyRTFilter.argtypes = [pRvRTFilter]
rvDestroyRTFilter.restype = None

rvRTFilterDelay = dll.rvRTFilterDelay
rvRTFilterDelay.argtypes = [ctypes.c_int]
rvRTFilterDelay.restype = ctypes.c_int

rvRTFilterDelayed = dll.rvRTFilterDelayed
rvRTFilterDelayed.argtypes = [pRvRTFilter]
rvRTFilterDelayed.restype = ctypes.c_int

rvRTFilterMaxOutputSize = dll.rvRTFilterMaxOutputSize
rvRTFilterMaxOutputSize.argtypes = [ctypes.c_int, ctypes.c_int]
rvRTFilterMaxOutputSize.restype = ctypes.c_int

class Procressor:
    def __init__(self, kernel, maxNX):
        self.kernel = kernel

        kernelSize = len(self.kernel)
        self.delay = kernelSize // 2
        if(kernelSize < 3):
            raise ValueError("length of kernel cannot be less than 3")
        if(kernelSize % 2 == 0):
            raise ValueError("length of kernel must be odd")
        
        self.proc = rvCreateRTFilter(kernel, kernelSize, maxNX)
    
    def __del__(self):
        rvDestroyRTFilter(self.proc)

    @property
    def delayed(self):
        return rvRTFilterDelayed(self.proc)

    def __call__(self, x):
        nX = 0 if(x is None) else len(x)
        nOut = rvRTFilterNextOutputSize(self.proc, nX)
        out = np.zeros(nOut, dtype = np.float64)
        if(nX == 0):
            v = rvCallRTFilter.argtypes
            rvCallRTFilter.argtypes = [pRvRTFilter, ctypes.POINTER(RvReal), ctypes.c_int, RvReal_1d]
            realOut = rvCallRTFilter(self.proc, None, 0, out)
            rvCallRTFilter.argtypes = v
        else:
            realOut = rvCallRTFilter(self.proc, x, nX, out)
        if(nOut == 0):
            return None
        else:
            return out