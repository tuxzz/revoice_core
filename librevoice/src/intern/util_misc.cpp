#define _USE_MATH_DEFINES
#include "util_p.hpp"

#include "find.hpp"
#include <cmath>

using namespace ReVoice;

void rvExitCheck()
{ rvCheckAllocated(); }

FrameRange rvGetFrameRange(int inputLen, int center, int size)
{
  rvAssert(inputLen > 0, "inputLen must be greater than 0");
  rvAssert(size > 0, "size must be greater than 0");
  int leftSize = size / 2;
  int rightSize = size - leftSize;

  int inputBegin = std::min(inputLen, std::max(center - leftSize, 0));
  int inputEnd = std::max(0, std::min(center + rightSize, inputLen));
  int outBegin = std::max(leftSize - center, 0);
  int outEnd = outBegin + (inputEnd - inputBegin);

  return {inputBegin, inputEnd, outBegin, outEnd};
}

void rvGetFrame(const RvReal *x, int nX, int center, int size, RvReal *out)
{
  rvAssert(x, "x cannot be nullptr");
  rvAssert(out, "out cannot be nullptr");
  rvAssert(!isInPlace(x, nX, out), "in-place operation is unsupported.");
  auto fr = rvGetFrameRange(nX, center, size);
  std::fill(out, out + size, 0.0);
  std::copy(x + fr.ib, x + fr.ie, out + fr.ob);
}

int rvGetNFrame(int inputSize, int hopSize)
{
  return inputSize % hopSize != 0 ? inputSize / hopSize + 1 : inputSize / hopSize;
}

void rvTailPad(const RvReal *x, int nX, RvReal *out, int nOut)
{
  rvAssert(x && out, "x or out cannot be nullptr");
  rvAssert(nX > 0 && nX <= nOut, "nX must be in range [1, nOut]");
  rvAssert(!isInPlace(x, nX, out), "in-place operation is unsupported.");
  std::copy(x, x + nX, out);
  std::fill(out + nX, out + nOut, 0.0);
}

void rvHeadPad(const RvReal *x, int nX, RvReal *out, int nOut)
{
  rvAssert(x && out, "x or out cannot be nullptr");
  rvAssert(nX > 0 && nX <= nOut, "nX must be in range [1, nOut]");
  rvAssert(!isInPlace(x, nX, out), "in-place operation is unsupported.");
  std::fill(out, out + nOut - nX, 0.0);
  std::copy(x, x + nX, out + nOut - nX);
}

void rvZeroPad(const RvReal *x, int nX, RvReal *out, int nOut)
{
  rvAssert(x && out, "x or out cannot be nullptr");
  rvAssert(nX > 0 && nX <= nOut, "nX must be in range [1, nOut]");
  rvAssert(nX % 2 == 0, "nX must be even");
  rvAssert(!isInPlace(x, nX, out), "in-place operation is unsupported.");

  int halfSize = nX / 2;
  std::copy(x, x + halfSize, out + nOut - halfSize);
  std::copy(x + halfSize, x + nX, out);
  std::fill(out + halfSize, out + nOut - halfSize, 0.0);
}

void rvFirwinSingleBand(int numtaps, RvReal passBegin, RvReal passEnd, const char *window, bool scale, RvReal nyq, RvReal *out)
{
  rvAssert(numtaps > 0, "numtaps must be greater than 0");
  rvAssert(passBegin >= 0.0 && passBegin <= passEnd, "passBegin must be in range [0.0, passEnd]");
  rvAssert(window, "window must be specified");
  rvAssert(nyq > 0.0, "nyq must be greater than 0");
  rvAssert(passEnd <= nyq, "passEnd or passBegin cannot be greater than nyq");
  rvAssert(out, "out cannot be nullptr");
  rvAssert(std::abs(passEnd - nyq) > 1e-32, "A filter with an even number of coefficients must have zero response at the Nyquist rate");

  auto windowInfo = rvGetWindow(window);
  rvAssert(windowInfo, "unsupported window");

  passBegin /= nyq;
  passEnd /= nyq;

  windowInfo->func(numtaps, out);
  RvReal alpha = 0.5 * static_cast<RvReal>(numtaps - 1);
  for(int i = 0; i < numtaps; ++i)
  {
    RvReal m = static_cast<RvReal>(i) - alpha;
    out[i] *= passEnd * sinc(passEnd * m) - passBegin * sinc(passBegin * m);
  }

  if(scale)
  {
    double f;
    if(passBegin == 0.0)
      f = 0.0;
    else if(passEnd == 1.0)
      f = 1.0;
    else
      f = 0.5 * (passBegin, passEnd);
    RvReal scaleFac = 0.0;
    f *= M_PI;
    for(int i = 0; i < numtaps; ++i)
    {
      RvReal m = static_cast<RvReal>(i) - alpha;
      scaleFac += std::cos(f * m) * out[i];
    }
    for(int i = 0; i < numtaps; ++i)
      out[i] /= scaleFac;
  }
}

ParabolicInterpolationResult rvParabolicInterp(const RvReal *x, int nX, int i, bool overAdjust)
{
  rvAssert(x, "x cannot be nullptr");
  rvAssert(nX, "nX must be greater than 0");
  rvAssert(i >= 0 && i < nX, "i must be in range [0, nX)");

  if(i > 0 && i < nX - 1)
  {
    RvReal s0 = x[i - 1];
    RvReal s1 = x[i];
    RvReal s2 = x[i + 1];
    RvReal a = (s0 + s2) / 2.0 - s1;
    if(a < 1e-32)
      return {static_cast<RvReal>(i), x[i]};
    RvReal b = s2 - s1 - a;
    RvReal adjustment = -(b / a * 0.5);
    if(!overAdjust && std::abs(adjustment) > 1.0)
      adjustment = 0.0;
    RvReal x = static_cast<RvReal>(i) + adjustment;
    RvReal y = a * adjustment * adjustment + b * adjustment + s1;
    return {x, y};
  }
  return {static_cast<RvReal>(i), x[i]};
}

void rvDoLinearInterp1d(const RvReal *xSrc, const RvReal *ySrc, int yRowSize, int nSrc, const RvReal *xIpl, int nIpl, RvReal *out)
{
  rvAssert(xSrc, "xSrc cannot be nullptr");
  rvAssert(ySrc, "ySrc cannot be nullptr");
  rvAssert(xIpl, "xIpl cannot be nullptr");
  rvAssert(out, "out cannot be nullptr");
  rvAssert(yRowSize > 0, "yRowSize must be greater than 0");
  rvAssert(nSrc > 0, "nSrc must be greater than 0");
  rvAssert(nIpl > 0, "nIpl must be greater than 0");

  int cachedIdx = 0;
  for(int iIpl = 0; iIpl < nIpl; ++iIpl)
  {
    int iSrc;
    if(xSrc[cachedIdx] <= xIpl[iIpl])
      iSrc = cachedIdx;
    else
      iSrc = static_cast<int>(matchLEQ(xSrc, xSrc + nSrc, xIpl[iIpl]) - xSrc);

    if(iSrc == nSrc - 1)
      out[iIpl] = xSrc[iSrc];
    else if(iIpl <= xSrc[0])
      out[iIpl] = xSrc[0];
    else
      out[iIpl] = lerp(xSrc[iSrc], xSrc[iSrc + 1], (xIpl[iIpl] - xSrc[iSrc]) / (xSrc[iSrc + 1] - xSrc[iSrc]));
  }
}

void rvSimpleDCRemove(RvReal *xo, int n)
{
  rvAssert(xo || n == 0, "xo must be greater than 0");
  rvAssert(n >= 0, "n must cannot be less than 0");
  RvReal m = mean(xo, n);
  for(int i = 0; i < n; ++i)
    xo[i] -= m;
}