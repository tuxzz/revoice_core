import sys
import numpy as np

from . import hmm

class SparseHMM:
    def __init__(self, init, frm, to, transProb, nMaxBackward):
        self.proc = hmm.SparseHMM(init, frm, to, transProb)

        self.init = init
        self.frm = frm
        self.to = to
        self.transProb = transProb

        self.nMaxBackward = int(nMaxBackward)
        assert nMaxBackward > 0

        nState = len(self.init)
        self.oldDelta = None
        self.psi = []
    
    def feed(self, obs):
        nState = len(self.init)

        if(self.oldDelta is None):
            self.oldDelta = self.init * obs
            deltaSum = np.sum(self.oldDelta)
            if(deltaSum > 0.0):
                self.oldDelta /= deltaSum
            self.psi.append(np.zeros(nState, dtype = np.int))
        else:
            psiFrame = np.zeros(nState, dtype = np.int)
            delta, psiFrame = self.proc.viterbiForwardRest(obs, self.oldDelta)
            deltaSum = np.sum(delta)
            self.psi.append(psiFrame)
            if(len(self.psi) > self.nMaxBackward):
                self.psi = self.psi[-self.nMaxBackward:]

            if(deltaSum > 0.0):
                self.oldDelta = delta / deltaSum
            else:
                print("WARNING: Viterbi decoder has been fed some zero probabilities.", file = sys.stderr)
                self.oldDelta.fill(1.0 / nState)
    
    def viterbiDecode(self, nBackward):
        # init backward step
        bestStateIdx = np.argmax(self.oldDelta)
        path = np.ndarray(nBackward, dtype = np.int) # the final output path
        path[-1] = bestStateIdx

        # rest of backward step
        localPsi = self.psi[-nBackward:]
        for iFrame in reversed(range(nBackward - 1)):
            path[iFrame] = localPsi[iFrame + 1][path[iFrame + 1]]
        return path