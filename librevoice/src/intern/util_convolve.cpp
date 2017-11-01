#include "util_p.hpp"

using namespace ReVoice;

struct RvFFTConvolver
{
  RvRFFT *rfft;
  RvIRFFT *irfft;
  RvComplex *cWorkMem;
  RvReal *rWorkMem;
  int maxSize;
};

void rvConvolve(const RvReal *x, int nX, const RvReal *y, int nY, RvReal *out)
{
  rvAssert(x, "x cannot be nullptr");
  rvAssert(y, "y cannot be nullptr");
  rvAssert(out, "out cannot be nullptr");
  rvAssert(nX > 0 && nY > 0, "nX or nY must be greater than 0");
  for(int i = 0; i < nX + nY - 1; ++i)
  {
    RvReal s = 0.0;
    int jBegin = std::max(0, i - nX + 1);
    int jEnd = std::min(i + 1, nY);
    for(int j = jBegin; j < jEnd; ++j)
      s += y[j] * x[i - j];
    out[i] = s;
  }
}

RvFFTConvolver *rvCreateFFTConvolver(int maxSize)
{
  rvAssert(maxSize > 0, "maxSize must be greater than 0");
  rvAssert(roundUpToPowerOf2(maxSize) == maxSize, "maxSize must be power of 2");
  int nF = maxSize / 2 + 1;
  auto convolver = new RvFFTConvolver;
  convolver->rfft = rvCreateRFFT(maxSize);
  convolver->irfft = rvCreateIRFFT(maxSize);
  convolver->cWorkMem = RVALLOC(RvComplex, nF * 2);
  convolver->rWorkMem = RVALLOC(RvReal, maxSize);
  convolver->maxSize = maxSize;
  return convolver;
}

void rvDestroyFFTConvolver(RvFFTConvolver *convolver)
{
  rvAssert(convolver, "convolver cannot be nullptr");
  rvFree(convolver->rWorkMem);
  rvFree(convolver->cWorkMem);
  rvDestroyRFFT(convolver->rfft);
  rvDestroyIRFFT(convolver->irfft);
  delete convolver;
}

void rvFFTConvolve(RvFFTConvolver *convolver, const RvReal *x, int nX, const RvReal *y, int nY, RvReal *out)
{
  int nOut = nX + nY - 1;
  int nF = convolver->maxSize / 2 + 1;
  rvAssert(convolver, "convolver cannot be nullptr");
  rvAssert(x, "x cannot be nullptr");
  rvAssert(y, "y cannot be nullptr");
  rvAssert(out, "out cannot be nullptr");
  rvAssert(nX > 0 && nY > 0, "nX or nY must be greater than 0");
  rvAssert(nOut <= convolver->maxSize, "nX + nY - 1 cannot be greater than maxSize");

  auto a = convolver->cWorkMem;
  auto b = convolver->cWorkMem + nF;
  rvTailPad(x, nX, convolver->rWorkMem, convolver->maxSize);
  rvDoRFFT(convolver->rfft, convolver->rWorkMem, a);
  rvTailPad(y, nY, convolver->rWorkMem, convolver->maxSize);
  rvDoRFFT(convolver->rfft, convolver->rWorkMem, b);
  for(int i = 0; i < nF; ++i)
    a[i] *= b[i];
  rvDoIRFFT(convolver->irfft, a, convolver->rWorkMem);
  std::copy(convolver->rWorkMem, convolver->rWorkMem + nOut, out);
}

int rvGetFFTConvolverSize(const RvFFTConvolver *convolver)
{
  rvAssert(convolver, "convolver cannot be nullptr");
  return convolver->maxSize;
}