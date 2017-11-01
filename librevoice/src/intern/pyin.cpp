#include "../pyin.h"

#include "util_p.hpp"
#include <cmath>

using namespace ReVoice;

void rvPYinNormalizedPdf(RvReal a, RvReal b, RvReal begin, RvReal end, RvReal *out, int n)
{
  RvReal step = (end - begin) / static_cast<RvReal>(n);
  for(int x = 0; x < n; ++x)
  {
    RvReal rx = static_cast<RvReal>(x) * step;
    out[x] = std::pow(rx, a - 1.0) * std::pow(1.0 - rx, b - 1.0);
  }

  for(int i = n - 2; i >= 0; --i)
  {
    if(out[i] < out[i + 1])
      out[i] = out[i + 1];
  }

  RvReal s = sum(out, n);
  for(int i = 0; i < n; ++i)
    out[i] /= s;
}