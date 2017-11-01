#pragma once

#include "util.h"

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct RvYinDifferenceWorker RvYinDifferenceWorker;

typedef struct RvYinProcessorParameter
{
  RvReal samprate;
  RvReal minFreq, maxFreq;
  RvReal valleyThreshold, valleyStep;
  int hopSize, windowSize;
  bool prefilter;
} RvYinProcessorParameter;

RV_EXPORT RvYinDifferenceWorker *rvCreateYinDifferenceWorker(int maxNX);
RV_EXPORT void rvDestroyYinDifferenceWorker(RvYinDifferenceWorker *worker);
RV_EXPORT void rvYinDoDifference(RvYinDifferenceWorker *worker, const RvReal *x, int nX, RvReal *out);

RV_EXPORT void rvYinCumulativeDifference(RvReal *xo, int n);
RV_EXPORT int rvYinFindValleys(const RvReal *x, int nX, RvReal minFreq, RvReal maxFreq, RvReal sr, RvReal threshold, RvReal step, int *out, int maxNOut);
RV_EXPORT void rvYinDoPrefilter(RvReal *xo, int n, RvReal maxFreq, RvReal sr);

RV_EXPORT RvYinProcessorParameter *rvCreateYinProcessorParameter(RvReal minFreq, RvReal maxFreq, RvReal sr);
RV_EXPORT void rvDestroyYinProcessorParameter(RvYinProcessorParameter *param);
RV_EXPORT void rvCallYin(const RvYinProcessorParameter *param, const RvReal *x, int nX, bool removeDC, RvReal *out);

#ifdef __cplusplus
}
#endif