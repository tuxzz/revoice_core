#pragma once

#include "../hmm.h"

typedef struct RvSparseHMM
{
  RvReal *init;
  int *frm, *to;
  RvReal *transProb;
  int nState, nTrans;

  RvReal *transTemp;
} RvSparseHMM;

namespace ReVoice
{
  void ctorSparseHMM(RvSparseHMM &self, RvReal *init, int *frm, int *to, RvReal *transProb, int nState, int nTrans);
  void dtorSparseHMM(RvSparseHMM &self);
} // namespace ReVoice