#pragma once

#include "util.h"

#ifdef __cplusplus
extern "C"
{
#endif

RV_EXPORT void rvPYinNormalizedPdf(RvReal a, RvReal b, RvReal begin, RvReal end, RvReal *out, int n);

#ifdef __cplusplus
}
#endif