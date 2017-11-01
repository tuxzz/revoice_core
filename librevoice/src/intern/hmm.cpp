#include "hmm_p.hpp"

#include "./util_p.hpp"

void rvHMMViterbiForwardRest(const RvSparseHMM *self, RvReal *oldDelta, RvReal *obs, RvReal *newDelta, int *psi)
{
  auto fromState = self->frm;
  auto toState = self->to;
  auto transProb = self->transProb;
  int nState = self->nState;
  int nTrans = self->nTrans;

  auto currValue = self->transTemp;
  for(int i = 0; i < nTrans; ++i)
    currValue[i] = oldDelta[fromState[i]] * transProb[i];
  std::fill(newDelta, newDelta + nState, 0.0);
  std::fill(psi, psi + nState, 0);
  for(int iTrans = 0; iTrans < nTrans; ++iTrans)
  {
    int ts = toState[iTrans];
    if(currValue[iTrans] > newDelta[ts])
    {
      newDelta[ts] = currValue[iTrans];
      psi[ts] = fromState[iTrans];
    }
  }

  for(int iState = 0; iState < nState; ++iState)
    newDelta[iState] *= obs[iState];
}

namespace ReVoice
{
  void ctorSparseHMM(RvSparseHMM &self, RvReal *init, int *frm, int *to, RvReal *transProb, int nState, int nTrans)
  {
    self.init = RVALLOC(RvReal, nState);
    self.frm = RVALLOC(int, nTrans);
    self.to = RVALLOC(int, nTrans);
    self.transProb = RVALLOC(RvReal, nTrans);
    self.transTemp = RVALLOC(RvReal, nTrans);

    self.nState = nState;
    self.nTrans = nTrans;

    std::copy(init, init + nState, self.init);
    std::copy(frm, frm + nTrans, self.frm);
    std::copy(to, to + nTrans, self.to);
    std::copy(transProb, transProb + nTrans, self.transProb);
  }

  void dtorSparseHMM(RvSparseHMM &self)
  {
    rvFree(self.init);
    rvFree(self.frm);
    rvFree(self.to);
    rvFree(self.transProb);
    rvFree(self.transTemp);
  }
} // namespace ReVoice