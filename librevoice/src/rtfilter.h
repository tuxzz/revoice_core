#pragma once
#include "util.h"

#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus

typedef struct RvRTFilter RvRTFilter;

RV_EXPORT RvRTFilter *rvCreateRTFilter(const RvReal *kernel, int kernelSize, int maxNX);
RV_EXPORT int rvRTFilterNextOutputSize(const RvRTFilter *rtfilter, int nX);
RV_EXPORT int rvCallRTFilter(RvRTFilter *rtfilter, const RvReal *x, int nX, RvReal *out);
RV_EXPORT void rvDestroyRTFilter(RvRTFilter *rtfilter);

RV_EXPORT int rvRTFilterDelay(int kernelSize);
RV_EXPORT int rvRTFilterDelayed(const RvRTFilter *rtfilter);
RV_EXPORT int rvRTFilterMaxOutputSize(int maxNX, int kernelSize);

#ifdef __cplusplus
}
#endif // __cplusplus
