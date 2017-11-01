#pragma once

#include "util.h"

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct RvSparseHMM RvSparseHMM;
RV_EXPORT void rvHMMViterbiForwardRest(const RvSparseHMM *sparseHMM, RvReal *oldDelta, RvReal *obs, RvReal *newDelta, int *psi);

#ifdef __cplusplus
}
#endif