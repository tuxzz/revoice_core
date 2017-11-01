#include "../rtpyin.h"

#include "util_p.hpp"
#include "../yin.h"
#include "../pyin.h"
#include "../rtfilter.h"

using namespace ReVoice;

typedef struct RvRTPYinProcessor
{
  RvRTPYinProcessorParameter param;
  RvRTFilter *filterProc;
  RvYinDifferenceWorker *differenceWorker;
  RvReal *buffer, *differenceTemp;
  int internalDelayed, bufferUsed;
  int bufferSize;
} RvRTPYinProcessor;

RvRTPYinProcessorParameter *rvCreateRTPYinProcessorParameter(RvReal minFreq, RvReal maxFreq, RvReal sr, RvReal *pdf, int pdfSize)
{
  RvReal nyq = sr / 2.0;
  rvAssert(sr > 0.0, "invalid samprate");
  rvAssert(maxFreq < nyq && maxFreq > minFreq && minFreq > 0.0, "invalid minFreq or maxFreq");
  rvAssert(pdf || pdfSize == 0, "pdf cannot be nullptr with non-zero pdfSize");

  bool selfAllocPdf = !pdf;
  if(selfAllocPdf)
  {
    pdf = RVALLOC(RvReal, 128);
    pdfSize = 128;
    rvPYinNormalizedPdf(1.7, 6.8, 0.0, 1.0, pdf, 128);
  }

  int hopSize = static_cast<int>(roundUpToPowerOf2(sr * 0.0025));

  return new RvRTPYinProcessorParameter({
    pdf,
    sr,
    minFreq, maxFreq,
    1.0, 0.01,
    0.02, 5.0,
    1.0,
    hopSize, std::max(static_cast<int>(roundUpToPowerOf2(sr / minFreq * 4.0)), hopSize),
    4, pdfSize,
    hopSize,
    true, selfAllocPdf
  });
}

void rvDestroyRTPYinProcessorParameter(RvRTPYinProcessorParameter *param)
{
  if(param->isPdfDefault)
    rvFree(param->pdf);
  delete param;
}

const RvRTPYinProcessorParameter *rvRTPYinParam(const RvRTPYinProcessor *self)
{ return &(self->param); }

RvRTPYinProcessor *rvCreateRTPYinProcessor(const RvRTPYinProcessorParameter *param)
{
  RvReal nyq = param->samprate / 2.0;
  rvAssert(param->pdf && param->pdfSize > 0, "pdf cannot be nullptr and pdfSize must be greater than 0");
  rvAssert(param->samprate > 0.0, "samprate must be greater than 0");
  rvAssert(param->maxFreq < nyq && param->maxFreq > param->minFreq && param->minFreq > 0.0, "invalid minFreq or maxFreq");
  rvAssert(param->valleyThreshold >= 0.0 && param->valleyStep > 0.0, "invalid valleyThreshold or valleyStep");
  rvAssert(param->probThreshold >= 0.0 && param->weightPrior > 0.0, "invalid probThreshold or weightPrior");
  rvAssert(param->bias > 0.0, "invalid bias");
  rvAssert(param->hopSize > 0 && param->maxWindowSize >= param->hopSize, "invalid hopSize or maxWindowSize");
  rvAssert(param->maxIter >= 1 && param->pdfSize > 0, "invalid maxIter or pdfSize");
  rvAssert(param->maxInputSegment > 0, "invalid maxInputSegment");

  auto self = new RvRTPYinProcessor;
  self->param = *param;
  auto pdf = self->param.pdf;
  self->param.pdf = RVALLOC(RvReal, self->param.pdfSize);
  std::copy(pdf, pdf + self->param.pdfSize, self->param.pdf);

  int filterOrder = 1;
  if(param->prefilter)
  {
    filterOrder = static_cast<int>(2048.0 * param->samprate / 44100.0);
    if(filterOrder % 2 == 0)
      filterOrder += 1;
    auto kernel = RVALLOC(RvReal, filterOrder);
    rvFirwinSingleBand(filterOrder, 0.0, std::max(1500.0, param->maxFreq * 4.0), "blackman", true, param->samprate / 2.0, kernel);
    self->filterProc = rvCreateRTFilter(kernel, filterOrder, param->maxInputSegment);
    rvFree(kernel);
  }
  self->differenceWorker = rvCreateYinDifferenceWorker(param->maxWindowSize);
  self->internalDelayed = false;
  self->bufferUsed = self->param.maxWindowSize / 2;

  self->bufferSize = param->maxWindowSize + std::max(param->maxInputSegment, rvRTFilterMaxOutputSize(param->maxInputSegment, filterOrder));
  self->buffer = RVALLOC(RvReal, self->bufferSize);
  std::fill(self->buffer, self->buffer + self->param.maxWindowSize / 2, 0.0);
  std::fill(self->buffer + self->param.maxWindowSize / 2, self->buffer + self->bufferSize, 0.0);

  self->differenceTemp = RVALLOC(RvReal, param->maxWindowSize / 2);

  return self;
}

int rvRTPYinDelayed(const RvRTPYinProcessor *self)
{ return self->param.prefilter ? self->internalDelayed + rvRTFilterDelayed(self->filterProc) : self->internalDelayed; }

int rvRTPYinBufferUsed(RvRTPYinProcessor *self)
{ return self->bufferUsed; }

void rvRTPYinDumpBuffer(RvRTPYinProcessor *self, RvReal *out)
{ std::copy(self->buffer, self->buffer + self->bufferUsed, out); }

int rvCallRTPYin(RvRTPYinProcessor *self, const RvReal *x, int nX, RvReal *out, int maxOut)
{
  rvAssert(x || nX == 0, "x cannot be nullptr with non-zero nX");
  rvAssert(nX >= 0 && nX <= self->param.maxInputSegment, "invalid nX");
  rvAssert(maxOut > 0 && maxOut <= 128, "maxOut must be in range (0, 128]");
  rvAssert(out, "out cannot be nullptr");
  int halfMaxWindowSize = self->param.maxWindowSize / 2;

  int nAppended;
  if(self->param.prefilter)
    nAppended = rvCallRTFilter(self->filterProc, x, nX, self->buffer + self->bufferUsed);
  else
  {
    std::copy(x, x + nX, self->buffer + self->bufferUsed);
    nAppended = nX;
  }
  if(nAppended == 0)
  {
    if(self->internalDelayed > 0)
    {
      std::fill(self->buffer + self->bufferUsed, self->buffer + self->bufferUsed + self->param.hopSize, 0.0);
      self->bufferUsed += self->param.hopSize;
    }
    else
      return -1;
  }
  else
  {
    self->internalDelayed += nAppended;
    self->bufferUsed += nAppended;
    if(self->bufferUsed < self->param.maxWindowSize)
      return -1;
  }

  rvAssert(self->bufferUsed >= self->param.maxWindowSize, "internal error");

  /* do pyin */
  int windowSize = 0;
  int newWindowSize = std::max(static_cast<int>(roundUpToPowerOf2(self->param.samprate / self->param.minFreq * 4.0)), self->param.hopSize * 2);
  int iIter = 0;

  int valleyIndexList[127];
  int nValley;
  while(newWindowSize != windowSize && iIter < self->param.maxIter)
  {
    windowSize = newWindowSize;
    int halfDelta = (self->param.maxWindowSize - windowSize) / 2;
    RvReal *frame = self->buffer + halfDelta;
    rvYinDoDifference(self->differenceWorker, frame, windowSize, self->differenceTemp);
    rvYinCumulativeDifference(self->differenceTemp, windowSize / 2);
    nValley = rvYinFindValleys(self->differenceTemp, windowSize / 2, self->param.minFreq, self->param.maxFreq, self->param.samprate, self->param.valleyThreshold, self->param.valleyStep, valleyIndexList, 127);
    if(nValley > 0)
    {
      RvReal possibleFreq = clip(self->param.minFreq, self->param.samprate / valleyIndexList[nValley - 1] - 20.0, self->param.maxFreq);
      newWindowSize = std::max(static_cast<int>(std::ceil(self->param.samprate / possibleFreq * 4.0)), self->param.hopSize * 2);
      if(newWindowSize % 2 != 0)
        newWindowSize += 1;
      iIter += 1;
    }
  }

  RvReal probTotal = 0.0;
  RvReal weightedProbTotal = 0.0;
  for(int iValley = 0; iValley < nValley; ++iValley)
  {
    auto result = rvParabolicInterp(self->differenceTemp, windowSize / 2, valleyIndexList[iValley], false);
    RvReal freq = self->param.samprate / result.x;
    RvReal v0 = iValley == 0 ? 1.0 : std::min(1.0, self->differenceTemp[valleyIndexList[iValley - 1]] + 1e-10);
    RvReal v1 = iValley == nValley - 1 ? 0.0 : std::max(0.0, self->differenceTemp[valleyIndexList[iValley + 1]]) + 1e-10;
    RvReal prob = 0.0;
    for(int i = static_cast<int>(v1 * self->param.pdfSize); i < static_cast<int>(v0 * self->param.pdfSize); ++i)
      prob += self->param.pdf[i] * (result.y < static_cast<RvReal>(i) / static_cast<RvReal>(self->param.pdfSize) ? 1.0 : 0.01);
    prob = std::min(prob, 0.99);
    prob *= self->param.bias;
    probTotal += prob;
    if(result.y < self->param.probThreshold)
      prob *= self->param.weightPrior;
    weightedProbTotal += prob;
    out[iValley * 2] = freq;
    out[iValley * 2 + 1] = prob;
  }

  if(nValley > 0 && weightedProbTotal != 0.0)
  {
    for(int iValley = 0; iValley < nValley; ++iValley)
      out[iValley * 2 + 1] *= probTotal / weightedProbTotal;
  }

  std::copy(self->buffer + self->param.hopSize, self->buffer + self->bufferUsed, self->buffer);
  std::fill(self->buffer + self->bufferUsed, self->buffer + self->bufferSize, 0.0);
  self->internalDelayed -= self->param.hopSize;
  self->bufferUsed -= self->param.hopSize;

  return nValley;
}

void rvDestroyRTPYinProcessor(RvRTPYinProcessor *self)
{
  rvFree(self->param.pdf);
  if(self->param.prefilter)
    rvDestroyRTFilter(self->filterProc);
  rvDestroyYinDifferenceWorker(self->differenceWorker);
  rvFree(self->buffer);
  rvFree(self->differenceTemp);
}

int rvRTPYinDelay(const RvRTPYinProcessorParameter *param)
{
  int delay = param->maxWindowSize / 2;
  if(param->prefilter)
  {
    int filterOrder = static_cast<int>(2048.0 * param->samprate / 44100.0);
    if(filterOrder % 2 == 0)
      filterOrder += 1;
    delay += filterOrder;
  }
  return delay;
}