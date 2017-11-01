#pragma once

#include "util.h"

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct RvRTPYinProcessorParameter RvRTPYinProcessorParameter;

typedef struct RvRTMonoPitchProcessorParameter
{
  RvReal samprate;
  RvReal maxTransSemitone, minFreq;
  RvReal transSelf, yinTrust, energyThreshold;
  int hopSize, nSemitone;
  int binPerSemitone, maxObsLength;
} RvRTMonoPitchProcessorParameter;
typedef struct RvRTMonoPitchProcessor RvRTMonoPitchProcessor;

RV_EXPORT RvRTMonoPitchProcessorParameter *rvCreateRTMonoPitchProcessorParameter(int hopSize, RvReal samprate, int nSemitone, RvReal maxTransSemitone, RvReal minFreq);
RV_EXPORT RvRTMonoPitchProcessorParameter *rvCreateRTMonoPitchProcessorParameterFromRTPYin(const RvRTPYinProcessorParameter *param);
RV_EXPORT void rvDestroyRTMonoPitchProcessorParameter(RvRTMonoPitchProcessorParameter *param);

RV_EXPORT RvRTMonoPitchProcessor *rvCreateRTMonoPitchProcessor(const RvRTMonoPitchProcessorParameter *param);
RV_EXPORT const RvRTMonoPitchProcessorParameter *rvRTMonoPitchParam(const RvRTMonoPitchProcessor *self);
RV_EXPORT void rvMonoPitchDumpObsTemp(const RvRTMonoPitchProcessor *self, RvReal *out);
RV_EXPORT int rvMonoPitchNextOutputLength(const RvRTMonoPitchProcessor *self);
RV_EXPORT int rvCallRTMonoPitch(RvRTMonoPitchProcessor *self, const RvReal *x, const RvReal *obsProb, int nObsProb, RvReal *out);
RV_EXPORT void rvDestroyRTMonoPitchProcessor(RvRTMonoPitchProcessor *self);

#ifdef __cplusplus
}
#endif