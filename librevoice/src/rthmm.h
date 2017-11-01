#pragma once

#include "util.h"

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct RvRTSparseHMM RvRTSparseHMM;

RV_EXPORT RvRTSparseHMM *rvCreateRTSparseHMM(RvReal *init, int *frm, int *to, RvReal *transProb, int nState, int nTrans, int nMaxBackward);
RV_EXPORT bool rvRTSparseHMMFeed(RvRTSparseHMM *rtSparseHMM, RvReal *obs);
RV_EXPORT int rvRTSparseHMMViterbiDecode(RvRTSparseHMM *rtSparseHMM, int *out, int nBackward);
RV_EXPORT int rvRTSparseHMMCurrentAvailable(RvRTSparseHMM *rtSparseHMM);
RV_EXPORT void rvDestroyRTSparseHMM(RvRTSparseHMM *rtSparseHMM);

#ifdef __cplusplus
}
#endif