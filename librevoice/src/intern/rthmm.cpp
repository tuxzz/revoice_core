#include "../rthmm.h"

#include "util_p.hpp"
#include "hmm_p.hpp"

#include <vector>

using namespace ReVoice;

typedef struct RvRTSparseHMM
{
  RvSparseHMM sparseHMM;

  RvReal *oldDelta, *deltaTemp;
  int *psiList;
  int nMaxBackward, psiUsed;
} RvRTSparseHMM;

RvRTSparseHMM *rvCreateRTSparseHMM(RvReal *init, int *frm, int *to, RvReal *transProb, int nState, int nTrans, int nMaxBackward)
{
  auto self = new RvRTSparseHMM;
  ctorSparseHMM(self->sparseHMM, init, frm, to, transProb, nState, nTrans);
  self->oldDelta = RVALLOC(RvReal, nState);
  self->deltaTemp = RVALLOC(RvReal, nState);
  self->psiList = RVALLOC(int, nMaxBackward * nState);
  self->nMaxBackward = nMaxBackward;
  self->psiUsed = 0;
  return self;
}

bool rvRTSparseHMMFeed(RvRTSparseHMM *self, RvReal *obs)
{
  rvAssert(obs, "obs cannot be nullptr");
  int nState = self->sparseHMM.nState;
  if(self->psiUsed == 0)
  {
    for(int i = 0; i < nState; ++i)
    {
      RvReal v = self->sparseHMM.init[i] * obs[i];
      self->oldDelta[i] = v;
    }
    RvReal deltaSum = sum(self->oldDelta, nState);
    if(deltaSum > 0.0)
    {
      for(int i = 0; i < nState; ++i)
        self->oldDelta[i] /= deltaSum;
    }
    std::fill(self->psiList + self->psiUsed * nState, self->psiList + (self->psiUsed + 1) * nState, 0);
    ++self->psiUsed;
    return true;
  }
  else
  {
    rvAssert(self->psiUsed <= self->nMaxBackward, "internal error");
    if(self->psiUsed == self->nMaxBackward)
    {
      std::copy(self->psiList + nState, self->psiList + self->psiUsed * nState, self->psiList);
      --self->psiUsed;
    }
    rvHMMViterbiForwardRest(&self->sparseHMM, self->oldDelta, obs, self->deltaTemp, self->psiList + self->psiUsed * nState);
    ++self->psiUsed;

    RvReal deltaSum = sum(self->deltaTemp, nState);
    if(deltaSum > 0.0)
    {
      for(int i = 0; i < nState; ++i)
        self->oldDelta[i] = self->deltaTemp[i] / deltaSum;
      return true;
    }
    else
    {
      warning("WARNING: Viterbi decoder has been fed some zero probabilities.");
      std::fill(self->oldDelta, self->oldDelta + nState, 1.0 / static_cast<RvReal>(nState));
      return false;
    }
  }
}

int rvRTSparseHMMViterbiDecode(RvRTSparseHMM *self, int *out, int nBackward)
{
  rvAssert(nBackward > 0 && nBackward <= self->nMaxBackward, "nBackward must be in range (0, nMaxBackward]");
  nBackward = std::min(self->psiUsed, nBackward);
  int nState = self->sparseHMM.nState;
  
  out[nBackward - 1] = argmax(self->oldDelta, nState);
  for(int i = nBackward - 2; i >= 0; --i)
     out[i] = self->psiList[(i + 1) * nState + out[i + 1]];

  return nBackward;
}

int rvRTSparseHMMCurrentAvailable(RvRTSparseHMM *self)
{ return self->psiUsed; }

void rvDestroyRTSparseHMM(RvRTSparseHMM *self)
{
  rvFree(self->psiList);
  rvFree(self->deltaTemp);
  rvFree(self->oldDelta);
  dtorSparseHMM(self->sparseHMM);
}