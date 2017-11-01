#include "../rtfilter.h"

#include "util_p.hpp"

using namespace ReVoice;

typedef struct RvRTFilter
{
  RvReal *kernel, *buffer;
  RvReal *convTemp;
  RvFFTConvolver *fftConv;
  int kernelSize, maxNX;
  int delayed;
} RvRTFilter;

RvRTFilter *rvCreateRTFilter(const RvReal *kernel, int kernelSize, int maxNX)
{
  rvAssert(kernelSize > 0, "kernelSize must be greater than 0");
  rvAssert(maxNX > 0, "maxNX must be greater than 0");
  rvAssert(kernel, "kernel cannot be nullptr");
  auto self = new RvRTFilter;
  self->kernel = RVALLOC(RvReal, kernelSize);
  self->buffer = RVALLOC(RvReal, kernelSize - 1);
  self->convTemp = RVALLOC(RvReal, kernelSize + maxNX - 1);
  self->fftConv = nullptr;
  self->kernelSize = kernelSize;
  self->maxNX = maxNX;
  self->delayed = 0;

  if(maxNX + kernelSize - 1 >= 128)
    self->fftConv = rvCreateFFTConvolver(roundUpToPowerOf2(maxNX + kernelSize - 1));

  std::copy(kernel, kernel + kernelSize, self->kernel);
  std::fill(self->buffer, self->buffer + kernelSize - 1, 0.0);

  return self;
}

int rvRTFilterNextOutputSize(const RvRTFilter *self, int nX)
{
  int delay = self->kernelSize / 2;
  if(nX == 0)
    return self->kernelSize / 2 - (delay - self->delayed);
  else if(self->delayed + nX <= delay)
    return 0;
  else
    return nX - (delay - (self->delayed + nX));
}

int rvCallRTFilter(RvRTFilter *self, const RvReal *x, int nX, RvReal *out)
{
  rvAssert(x || nX == 0, "x cannot be nullptr with non-zero nX");
  rvAssert(nX >= 0 && nX <= self->maxNX, "nX must be in range [0, maxNX]");
  int delay = self->kernelSize / 2;
  if(nX > 0)
  {
    int convedSize = nX + self->kernelSize - 1;
    if(nX + delay < 128)
      rvConvolve(x, nX, self->kernel, self->kernelSize, self->convTemp);
    else
      rvFFTConvolve(self->fftConv, x, nX, self->kernel, self->kernelSize, self->convTemp);
    for(int i = 0; i < self->kernelSize - 1; ++i)
      self->convTemp[i] += self->buffer[i];
    std::copy(self->convTemp + convedSize - self->kernelSize + 1, self->convTemp + convedSize, self->buffer);
    if(self->delayed < delay)
    {
      int nOut = 0;
      if(nX > delay - self->delayed)
      {
        std::copy(self->convTemp + delay - self->delayed, self->convTemp + nX, self->convTemp);
        nOut = nX - (delay - self->delayed);
      }
      self->delayed += std::min(nX, delay - self->delayed);
      std::copy(self->convTemp, self->convTemp + nOut, out);
      return nOut;
    }
    else
    {
      std::copy(self->convTemp, self->convTemp + nX, out);
      return nX;
    }
  }
  else
  {
    if(self->delayed == 0)
      return 0;
    std::copy(self->buffer + delay - self->delayed, self->buffer + self->kernelSize / 2, out);
    std::fill(self->buffer, self->buffer + self->kernelSize - 1, 0.0);
    int nOut = self->kernelSize / 2 - (delay - self->delayed);
    self->delayed = 0;
    return nOut;
  }  
}

void rvDestroyRTFilter(RvRTFilter *self)
{
  if(self->fftConv)
    rvDestroyFFTConvolver(self->fftConv);
  rvFree(self->convTemp);
  rvFree(self->buffer);
  rvFree(self->kernel);
}

int rvRTFilterDelay(int kernelSize)
{ return kernelSize - 1; }

int rvRTFilterDelayed(const RvRTFilter *self)
{ return self->delayed; }
int rvRTFilterMaxOutputSize(int maxNX, int kernelSize)
{ return std::max(maxNX, kernelSize / 2); }