#define _USE_MATH_DEFINES
#include "util_p.hpp"

#include <cmath>

using namespace ReVoice;

static RvWindowInfo windowList[] = {
  {"hanning", rvHanning, 1.5, 0.5},
  {"blackman", rvBlackman, 1.73, 0.42},
  {"blackmanharris", rvBlackmanHarris, 2.0044, (35875.0 - 3504.0 * M_PI) / 100000.0}
};

void rvHanning(int n, RvReal *out)
{
  rvAssert(n > 0, "n must be greater than 0");
  rvAssert(out, "out cannot be nullptr");
  RvReal a = 2.0 * M_PI / static_cast<RvReal>(n - 1);
  for(int i = 0; i < n; ++i)
    out[i] = 0.5 - 0.5 * std::cos(a * static_cast<RvReal>(i));
}

void rvBlackman(int n, RvReal *out)
{
  rvAssert(n > 0, "n must be greater than 0");
  rvAssert(out, "out cannot be nullptr");
  RvReal a = 2.0 * M_PI / static_cast<RvReal>(n - 1);
  RvReal b = 4.0 * M_PI / static_cast<RvReal>(n - 1);
  for(int i = 0; i < n; ++i)
  {
    auto ri = static_cast<RvReal>(i);
    out[i] = 0.42 - 0.5 * std::cos(a * ri) + 0.08 * std::cos(b * ri);
  }
}

void rvBlackmanHarris(int n, RvReal *out)
{
  rvAssert(n > 0, "n must be greater than 0");
  rvAssert(out, "out cannot be nullptr");
  RvReal a = 2.0 * M_PI / static_cast<RvReal>(n - 1);
  RvReal b = 4.0 * M_PI / static_cast<RvReal>(n - 1);
  RvReal c = 6.0 * M_PI / static_cast<RvReal>(n - 1);
  for(int i = 0; i < n; ++i)
  {
    auto ri = static_cast<RvReal>(i);
    out[i] = 0.35875 - 0.48829 * std::cos(a * ri) + 0.14128 * std::cos(b * ri) - 0.01168 * std::cos(c * ri);
  }
}

void rvGaussian(int n, RvReal stdev, RvReal *out)
{
  rvAssert(n > 0, "n must be greater than 0");
  rvAssert(!(std::isnan(stdev) || std::isinf(stdev)), "stdev cannot be nan or inf");
  rvAssert(out, "out cannot be nullptr");
  RvReal a = -0.5 * stdev * stdev;
  for(int i = 0; i < n; ++i)
    out[i] = std::exp(static_cast<RvReal>(i * i) / a);
}

int rvWindowCount()
{ return sizeof(windowList) / sizeof(RvWindowInfo); }

RvWindowInfo *rvWindowList()
{ return windowList; }

RvWindowInfo *rvGetWindow(const char *name)
{
  rvAssert(name, "name cannot be nullptr");
  for(int i = 0; i < rvWindowCount(); ++i)
  {
    if(!strcmp(name, windowList[i].name))
      return &(windowList[i]);
  }
  return nullptr;
}