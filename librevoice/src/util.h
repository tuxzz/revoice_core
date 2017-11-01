#pragma once
#include <stdbool.h>

typedef double RvReal;
#ifdef __cplusplus
#include <complex>
typedef std::complex<RvReal> RvComplex;
#else
#include <complex.h>
typedef RvReal complex RvComplex;
#endif // __cplusplus

#define RV_EXPORT __declspec(dllexport)

#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus

typedef struct RvRFFT RvRFFT;
typedef struct RvIRFFT RvIRFFT;
typedef struct RvFFTConvolver RvFFTConvolver;
typedef void (RvWindowFunc)(int n, RvReal *out);

typedef struct RvWindowInfo
{
  const char *name;
  RvWindowFunc *func;
  RvReal B, mean;
} RvWindowInfo;

typedef struct FrameRange
{
  int ib, ie, ob, oe;
} FrameRange;

typedef struct ParabolicInterpolationResult
{
  RvReal x, y;
} ParabolicInterpolationResult;

RV_EXPORT void rvExitCheck();

RV_EXPORT FrameRange rvGetFrameRange(int inputLen, int center, int size);
RV_EXPORT void rvGetFrame(const RvReal *x, int nX, int center, int size, RvReal *out);
RV_EXPORT int rvGetNFrame(int inputSize, int hopSize);
RV_EXPORT void rvTailPad(const RvReal *x, int nX, RvReal *out, int nOut);
RV_EXPORT void rvHeadPad(const RvReal *x, int nX, RvReal *out, int nOut);
RV_EXPORT void rvZeroPad(const RvReal *x, int nX, RvReal *out, int nOut);

RV_EXPORT RvRFFT *rvCreateRFFT(int fftSize);
RV_EXPORT void rvDestroyRFFT(RvRFFT *rfft);
RV_EXPORT void rvDoRFFT(RvRFFT *rfft, const RvReal *in, RvComplex *out);
RV_EXPORT int rvGetRFFTSize(const RvRFFT *rfft);

RV_EXPORT RvIRFFT *rvCreateIRFFT(int fftSize);
RV_EXPORT void rvDestroyIRFFT(RvIRFFT *irfft);
RV_EXPORT void rvDoIRFFT(RvIRFFT *irfft, const RvComplex *in, RvReal *out);
RV_EXPORT int rvGetIRFFTSize(const RvIRFFT *rfft);

RV_EXPORT void rvConvolve(const RvReal *x, int nX, const RvReal *y, int nY, RvReal *out);
RV_EXPORT RvFFTConvolver *rvCreateFFTConvolver(int maxSize);
RV_EXPORT void rvDestroyFFTConvolver(RvFFTConvolver *convolver);
RV_EXPORT void rvFFTConvolve(RvFFTConvolver *convolver, const RvReal *x, int nX, const RvReal *y, int nY, RvReal *out);
RV_EXPORT int rvGetFFTConvolverSize(const RvFFTConvolver *convolver);

RV_EXPORT void rvHanning(int n, RvReal *out);
RV_EXPORT void rvBlackman(int n, RvReal *out);
RV_EXPORT void rvBlackmanHarris(int n, RvReal *out);
RV_EXPORT void rvGaussian(int n, RvReal stdev, RvReal *out);

RV_EXPORT int rvWindowCount();
RV_EXPORT RvWindowInfo *rvWindowList();
RV_EXPORT RvWindowInfo *rvGetWindow(const char *name);

RV_EXPORT void rvFirwinSingleBand(int numtaps, RvReal passBegin, RvReal passEnd, const char *window, bool scale, RvReal nyq, RvReal *out);

RV_EXPORT ParabolicInterpolationResult rvParabolicInterp(const RvReal *x, int nX, int i, bool overAdjust);
RV_EXPORT void rvDoLinearInterp1d(const RvReal *xSrc, const RvReal *ySrc, int yRowSize, int nSrc, const RvReal *xIpl, int nIpl, RvReal *out);
RV_EXPORT void rvSimpleDCRemove(RvReal *xo, int n);

#ifdef __cplusplus
}
#endif // __cplusplus
