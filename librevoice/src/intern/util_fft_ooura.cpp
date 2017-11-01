#include "util_p.hpp"

#include "fftsg_h.h"

using namespace ReVoice;

struct RvRFFT
{
  bool inverse;
  int fftSize;
  double *transformBuffer;
};

RvRFFT *rvCreateRFFT(int fftSize)
{
  rvAssert(fftSize > 0, "fftSize must be greater than 0");
  rvAssert(roundUpToPowerOf2(fftSize) == fftSize, "fftSize must be power of 2");
  auto rfft = new RvRFFT;
  rfft->inverse = false;
  rfft->fftSize = fftSize;
  rfft->transformBuffer = RVALLOC(double, fftSize * 2);
  std::fill(rfft->transformBuffer, rfft->transformBuffer + fftSize * 2, 0.0);
  return rfft;
}

void rvDestroyRFFT(RvRFFT *rfft)
{
  rvAssert(rfft, "rfft cannot be nullptr");
  rvAssert(!rfft->inverse, "not a rfft object");
  rvFree(rfft->transformBuffer);
  delete rfft;
}

void rvDoRFFT(RvRFFT *rfft, const RvReal *in, RvComplex *out)
{
  rvAssert(rfft, "rfft cannot be nullptr");
  rvAssert(in, "in cannot be nullptr");
  rvAssert(out, "out cannot be nullptr");
  rvAssert(!rfft->inverse, "not a rfft object");

  for(int i = 0; i < rfft->fftSize; ++i)
  {
    rfft->transformBuffer[i * 2] = in[i];
    rfft->transformBuffer[i * 2 + 1] = 0.0;
  }
  cdft(rfft->fftSize * 2, -1, rfft->transformBuffer);
  for(int i = 0; i < rfft->fftSize / 2 + 1; ++i)
  {
    out[i].real(rfft->transformBuffer[i * 2]);
    out[i].imag(rfft->transformBuffer[i * 2 + 1]);
  }
}

int rvGetRFFTSize(const RvRFFT *rfft)
{
  rvAssert(rfft, "rfft cannot be nullptr");
  rvAssert(!rfft->inverse, "not a rfft object");
  return rfft->fftSize;
}

RvIRFFT *rvCreateIRFFT(int fftSize)
{
  rvAssert(fftSize % 2 == 0, "odd fftSize is not supported in IRFFT");
  rvAssert(fftSize > 0, "fftSize must be greater than 0");
  auto obj = rvCreateRFFT(fftSize);
  obj->inverse = true;
  return reinterpret_cast<RvIRFFT*>(obj);
}

void rvDestroyIRFFT(RvIRFFT *irfft)
{
  rvAssert(irfft, "irfft cannot be nullptr");
  rvAssert(reinterpret_cast<RvRFFT*>(irfft)->inverse, "not a irfft object");
  reinterpret_cast<RvRFFT*>(irfft)->inverse = false;
  rvDestroyRFFT(reinterpret_cast<RvRFFT*>(irfft));
}

void rvDoIRFFT(RvIRFFT *_irfft, const RvComplex *in, RvReal *out)
{
  rvAssert(_irfft, "irfft cannot be nullptr");
  rvAssert(in, "in cannot be nullptr");
  rvAssert(out, "out cannot be nullptr");
  rvAssert(reinterpret_cast<RvRFFT*>(_irfft)->inverse, "not a irfft object");

  RvRFFT *irfft = reinterpret_cast<RvRFFT*>(_irfft);

  for(int i = 0; i < irfft->fftSize / 2 + 1; ++i)
  {
    irfft->transformBuffer[i * 2] = in[i].real();
    irfft->transformBuffer[i * 2 + 1] = in[i].imag();
  }
  for(int i = irfft->fftSize / 2 + 1; i < irfft->fftSize; ++i)
  {
    irfft->transformBuffer[i * 2] = in[irfft->fftSize - i].real();
    irfft->transformBuffer[i * 2 + 1] = -in[irfft->fftSize - i].imag();
  }

  cdft(irfft->fftSize * 2, 1, irfft->transformBuffer);

  RvReal scaleFac = 1.0 / static_cast<double>(irfft->fftSize);
  for(int i = 0; i < irfft->fftSize; ++i)
    out[i] = irfft->transformBuffer[i * 2] * scaleFac;
}

int rvGetIRFFTSize(const RvIRFFT *irfft)
{
  rvAssert(irfft, "irfft cannot be nullptr");
  rvAssert(reinterpret_cast<const RvRFFT*>(irfft)->inverse, "not a irfft object");
  return reinterpret_cast<const RvRFFT*>(irfft)->fftSize;
}