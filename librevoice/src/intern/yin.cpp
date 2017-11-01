#include "../yin.h"

#include "util_p.hpp"

using namespace ReVoice;

typedef struct RvYinDifferenceWorker
{
  RvRFFT *rfft;
  RvIRFFT *irfft;
  RvComplex *cWorkMem;
  RvReal *rWorkMem;
  int nPadded;
} RvYinDifferenceWorker;

RvYinDifferenceWorker *rvCreateYinDifferenceWorker(int maxNX)
{
  rvAssert(maxNX > 0, "maxNX must be greater than 0");
  rvAssert(roundUpToPowerOf2(maxNX), "maxNX must be power of 2");
  int nTransformed = maxNX / 2 + 1;
  auto worker = new RvYinDifferenceWorker;
  worker->rfft = rvCreateRFFT(maxNX);
  worker->irfft = rvCreateIRFFT(maxNX);
  worker->cWorkMem = RVALLOC(RvComplex, 2 * nTransformed);
  worker->rWorkMem = RVALLOC(RvReal, maxNX);
  worker->nPadded = maxNX;
  return worker;
}

void rvDestroyYinDifferenceWorker(RvYinDifferenceWorker *worker)
{
  rvAssert(worker, "worker cannot be nullptr");
  rvFree(worker->rWorkMem);
  rvFree(worker->cWorkMem);
  rvDestroyIRFFT(worker->irfft);
  rvDestroyRFFT(worker->rfft);
  delete worker;
}

void rvYinDoDifference(RvYinDifferenceWorker *worker, const RvReal *x, int nX, RvReal *out)
{
  rvAssert(worker, "worker cannot be nullptr");
  rvAssert(x, "x cannot be nullptr");
  rvAssert(nX > 0, "nX must be greater than 0");
  rvAssert(out, "out cannot be nullptr");
  rvAssert(nX % 2 == 0, "nX must be even");

  int nPadded = worker->nPadded;
  int nOut = nX / 2;
  int nTransformed = nPadded / 2 + 1;

  rvAssert(nX <= nPadded, "nX cannot be greater than nPadded");

  // POWER TERM CALCULATION
  // ... for the power terms in equation (7) in the Yin paper
  out[0] = squareSum(x, nOut);
  for(int i = 1; i < nOut; ++i)
    out[i] = out[i - 1] - x[i - 1] * x[i - 1] + x[i + nOut] * x[i + nOut];

  // YIN-STYLE ACF via FFT
  auto buffer = worker->rWorkMem;
  auto transformedAudio = worker->cWorkMem;
  auto transformedKernel = worker->cWorkMem + nTransformed;

  rvTailPad(x, nX, buffer, nPadded);
  rvDoRFFT(worker->rfft, buffer, transformedAudio);

  for(int i = 0; i < nOut; ++i)
    buffer[i] = x[nOut - i - 1];
  std::fill(buffer + nOut, buffer + nPadded, 0.0);
  rvDoRFFT(worker->rfft, buffer, transformedKernel);

  for(int i = 0; i < nTransformed; ++i)
    transformedAudio[i] *= transformedKernel[i];
  rvDoIRFFT(worker->irfft, transformedAudio, buffer);

  auto powerTerm0 = out[0];
  for(int i = 0; i < nOut; ++i)
    out[i] += powerTerm0 - 2 * buffer[nOut - 1 + i];
}

void rvYinCumulativeDifference(RvReal *xo, int n)
{
  xo[0] = 1.0;
  RvReal s = 0.0;
  for(int i = 1; i < n; ++i)
  {
    s += xo[i];
    if(s == 0.0)
      xo[i] = 1.0;
    else
      xo[i] *= static_cast<RvReal>(i) / s;
  }
}

int rvYinFindValleys(const RvReal *x, int nX, RvReal minFreq, RvReal maxFreq, RvReal sr, RvReal threshold, RvReal step, int *out, int maxNOut)
{
  int iBegin = std::max(1, static_cast<int>(sr / maxFreq));
  int iEnd = std::min(nX - 1, static_cast<int>(std::ceil(sr / minFreq)));
  int nOut = 0;
  for(int i = iBegin; i < iEnd && nOut < maxNOut; ++i)
  {
    auto prev = x[i - 1];
    auto curr = x[i];
    auto next = x[i + 1];
    if(prev > curr && next > curr && curr < threshold)
    {
      threshold = curr - step;
      out[nOut] = i;
      ++nOut;
    }
  }

  return nOut;
}

void rvYinDoPrefilter(RvReal *xo, int n, RvReal maxFreq, RvReal sr)
{
  int filterOrder = static_cast<int>(2048.0 * sr / 44100.0);
  int halfFilterOrder = filterOrder / 2;
  if(filterOrder % 2 == 0)
    filterOrder += 1;
  auto buffer = RVALLOC(RvReal, n + filterOrder - 1);
  auto filter = RVALLOC(RvReal, filterOrder);
  rvFirwinSingleBand(filterOrder, 0.0, std::max(1250.0, maxFreq * 1.25), "blackman", true, sr / 2.0, filter);

  auto conv = rvCreateFFTConvolver(roundUpToPowerOf2(n + filterOrder - 1));
  rvFFTConvolve(conv, xo, n, filter, filterOrder, buffer);
  rvDestroyFFTConvolver(conv);
  rvFree(filter);

  for(int i = 0; i < n; ++i)
    xo[i] = buffer[halfFilterOrder + i];
  rvFree(buffer);
}

RvYinProcessorParameter *rvCreateYinProcessorParameter(RvReal minFreq, RvReal maxFreq, RvReal sr)
{
  rvAssert(sr > 0, "sr must be greater than 0");
  rvAssert(minFreq > 0, "minFreq must be greater than 0");
  rvAssert(maxFreq > minFreq, "maxFreq must be greater than minFreq");
  auto param = new RvYinProcessorParameter;
  param->samprate = sr;
  param->minFreq = minFreq;
  param->maxFreq = maxFreq;
  param->valleyThreshold = 0.5;
  param->valleyStep = 0.01;
  param->hopSize = roundUpToPowerOf2(static_cast<int>(std::round(sr * 0.0025)));
  param->windowSize = std::max(roundUpToPowerOf2(static_cast<int>(sr / minFreq) * 2), param->hopSize * 4);
  param->prefilter = true;
  return param;
}

void rvDestroyYinProcessorParameter(RvYinProcessorParameter *param)
{
  rvAssert(param, "param cannot be nullptr");
  delete param;
}

void rvCallYin(const RvYinProcessorParameter *param, const RvReal *x, int nX, bool removeDC, RvReal *out)
{
  rvAssert(param, "param cannot be nullptr");
  rvAssert(x, "x cannot be nullptr");
  rvAssert(nX > 0, "nX muse be greater than 0");
  rvAssert(out, "out cannot be nullptr");
  auto px = RVALLOC(RvReal, nX);
  std::copy(x, x + nX, px);
  if(removeDC)
    rvSimpleDCRemove(px, nX);
  if(param->prefilter)
    rvYinDoPrefilter(px, nX, param->maxFreq, param->samprate);

  int nHop = rvGetNFrame(nX, param->hopSize);
  auto frame = RVALLOC(RvReal, param->windowSize);
  auto buffer = RVALLOC(RvReal, param->windowSize / 2);
  auto worker = rvCreateYinDifferenceWorker(param->windowSize);
  int valleys[32];
  for(int iHop = 0; iHop < nHop; ++iHop)
  {
    rvGetFrame(px, nX, iHop * param->hopSize, param->windowSize, frame);
    rvYinDoDifference(worker, frame, param->windowSize, buffer);
    rvYinCumulativeDifference(buffer, param->windowSize / 2);
    int nValley = rvYinFindValleys(buffer, param->windowSize / 2, param->minFreq, param->maxFreq, param->samprate, param->valleyThreshold, param->valleyStep, valleys, 32);
    if(nValley > 0)
    {
      RvReal ipledX = rvParabolicInterp(buffer, param->windowSize / 2, valleys[nValley - 1], false).x;
      out[iHop] = param->samprate / ipledX;
    }
    else
      out[iHop] = 0.0;
  }
  rvDestroyYinDifferenceWorker(worker);
  rvFree(buffer);
  rvFree(frame);
  rvFree(px);
}