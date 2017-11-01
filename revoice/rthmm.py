import ctypes
import numpy as np
import numpy.ctypeslib as npct
from . import rtpyin

dll = ctypes.CDLL("librevoice.dll")
RvReal = ctypes.c_double
int_1d = npct.ndpointer(dtype = np.int, ndim = 1, flags = "C")
int_2d = npct.ndpointer(dtype = np.int, ndim = 2, flags = "C")
RvReal_1d = npct.ndpointer(dtype = np.float64, ndim = 1, flags = "C")
pRvReal = ctypes.POINTER(RvReal)
pint = ctypes.POINTER(ctypes.c_int)

class RvRTSparseHMM(ctypes.Structure):
    pass

pRvRTSparseHMM = ctypes.POINTER(RvRTSparseHMM)

rvCreateRTSparseHMM = dll.rvCreateRTSparseHMM
rvCreateRTSparseHMM.argtypes = [RvReal_1d, int_1d, int_1d, RvReal_1d, ctypes.c_int, ctypes.c_int, ctypes.c_int]
rvCreateRTSparseHMM.restype = pRvRTSparseHMM

rvRTSparseHMMPsiListUsed = dll.rvRTSparseHMMPsiListUsed
rvRTSparseHMMPsiListUsed.argtypes = [pRvRTSparseHMM]
rvRTSparseHMMPsiListUsed.restype = ctypes.c_int

rvRTSparseHMMDumpPsiList = dll.rvRTSparseHMMDumpPsiList
rvRTSparseHMMDumpPsiList.argtypes = [pRvRTSparseHMM, int_2d]
rvRTSparseHMMDumpPsiList.restype = ctypes.c_int

rvRTSparseHMMFeed = dll.rvRTSparseHMMFeed
rvRTSparseHMMFeed.argtypes = [pRvRTSparseHMM, RvReal_1d]
rvRTSparseHMMFeed.restype = ctypes.c_bool

rvRTSparseHMMViterbiDecode = dll.rvRTSparseHMMViterbiDecode
rvRTSparseHMMViterbiDecode.argtypes = [pRvRTSparseHMM, int_1d, ctypes.c_int]
rvRTSparseHMMViterbiDecode.restype = ctypes.c_int

rvRTSparseHMMCurrentAvailable = dll.rvRTSparseHMMCurrentAvailable
rvRTSparseHMMCurrentAvailable.argtypes = [pRvRTSparseHMM]
rvRTSparseHMMCurrentAvailable.restype = ctypes.c_int

rvDestroyRTSparseHMM = dll.rvDestroyRTSparseHMM
rvDestroyRTSparseHMM.argtypes = [pRvRTSparseHMM]
rvDestroyRTSparseHMM.restype = None

class SparseHMM:
    def __init__(self, init, frm, to, transProb, nMaxBackward):
        self.init = init
        self.frm = frm
        self.to = to
        self.transProb = transProb

        self.nMaxBackward = int(nMaxBackward)
        assert nMaxBackward > 0

        nState = len(self.init)
        nTrans = len(self.transProb)
        self.proc = rvCreateRTSparseHMM(init, frm, to, transProb, nState, nTrans, nMaxBackward)
    
    def __del__(self):
        rvDestroyRTSparseHMM(self.proc)

    @property
    def psi(self):
        n = rvRTSparseHMMPsiListUsed(self.proc)
        psi = np.zeros((n, self.proc.contents.sparseHMM.nState), dtype = np.int)
        rvRTSparseHMMDumpPsiList(self.proc, psi)
        return psi

    def feed(self, obs):
        rvRTSparseHMMFeed(self.proc, np.require(obs, np.float64, ("C_CONTIGUOUS",)))
    
    def viterbiDecode(self, nBackward):
        n = min(nBackward, rvRTSparseHMMCurrentAvailable(self.proc))
        out = np.zeros(n, dtype = np.int)
        realN = rvRTSparseHMMViterbiDecode(self.proc, out, nBackward)
        assert n == realN
        return out
    
    def _viterbiDecode(self, nBackward):
        # init backward step
        oldDelta = np.ctypeslib.as_array(self.proc.contents.oldDelta, shape = (self.proc.contents.sparseHMM.nState,))
        bestStateIdx = np.argmax(oldDelta)
        path = np.ndarray(nBackward, dtype = np.int) # the final output path
        path[-1] = bestStateIdx

        # rest of backward step
        localPsi = psi[-nBackward:]
        for iFrame in reversed(range(nBackward - 1)):
            path[iFrame] = localPsi[iFrame + 1][path[iFrame + 1]]
        return path