#pragma once

#include "util.h"

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct RvRTPYinProcessorParameter
{
  RvReal *pdf;
  RvReal samprate;
  RvReal minFreq, maxFreq;
  RvReal valleyThreshold, valleyStep;
  RvReal probThreshold, weightPrior;
  RvReal bias;
  int hopSize, maxWindowSize;
  int maxIter, pdfSize;
  bool prefilter, isPdfDefault;
} RvRTPYinProcessorParameter;

typedef struct RvRTPYinProcessor RvRTPYinProcessor;

RV_EXPORT RvRTPYinProcessorParameter *rvCreateRTPYinProcessorParameter(RvReal minFreq, RvReal maxFreq, RvReal sr, RvReal *pdf, int pdfSize);
RV_EXPORT void rvDestroyRTPYinProcessorParameter(RvRTPYinProcessorParameter *param);

RV_EXPORT RvRTPYinProcessor *rvCreateRTPYinProcessor(const RvRTPYinProcessorParameter *param);
RV_EXPORT const RvRTPYinProcessorParameter *rvRTPYinParam(const RvRTPYinProcessor *rtpyin);
RV_EXPORT int rvRTPYinDelayed(const RvRTPYinProcessor *rtpyin);
RV_EXPORT int rvCallRTPYin(RvRTPYinProcessor *rtpyin, const RvReal *x, int nX, RvReal *out, int maxOut);
RV_EXPORT int rvRTPYinBufferUsed(RvRTPYinProcessor *rtpyin);
RV_EXPORT void rvRTPYinDumpBuffer(RvRTPYinProcessor *rtpyin, RvReal *out);
RV_EXPORT void rvDestroyRTPYinProcessor(RvRTPYinProcessor *rtpyin);

RV_EXPORT int rvRTPYinDelay(const RvRTPYinProcessorParameter *param);

#ifdef __cplusplus
}
#endif