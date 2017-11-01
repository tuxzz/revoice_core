#pragma once
#include "../util.h"
#include <algorithm>
#include <numeric>

#include <cstdio>
#include <cstdlib>

#define _rv_S(a) _rv__S(a)
#define _rv__S(a) #a

#define rvAssert(a, ...) \
        (void)((!(a)) ?  ( \
            ( \
            fprintf(stderr, \
                "rvAssert failed: %s:%d, %s(), at \'%s\' " __VA_ARGS__ "\n", \
                __FILE__, __LINE__, __func__, _rv_S(a)), \
            std::abort(), \
            NULL)) : NULL)

#ifndef NDEBUG
#define rvDebugAssert(a, ...) \
        (void)((!(a)) ?  ( \
            ( \
            fprintf(stderr, \
                "rvDebugAssert failed: %s:%d, %s(), at \'%s\' " __VA_ARGS__ "\n", \
                __FILE__, __LINE__, __func__, _rv_S(a)), \
            std::abort(), \
            NULL)) : NULL)
#else
#define rvDebugAssert(a, ...)
#endif // NDEBUG

#define RVALLOC(t, n) ReVoice::rvAlloc< t >(n, __FILE__, __func__, __LINE__)

#if (defined(__AVX2__) || defined( __AVX__ ) || defined(_M_AMD64) || defined(_M_X64) || _M_IX86_FP == 2)
#ifndef __SSE2__
#define __SSE2__
#endif
#ifndef __SSE__
#define __SSE__
#endif
#endif

#if _M_IX86_FP == 1
#ifndef __SSE__
#define __SSE__
#endif
#endif

#if defined(_MSC_VER)
#define RV_ALIGNED(x) __declspec(align(x))
#elif defined(__GNUC__)
#define RV_ALIGNED(x) __attribute__ ((aligned(x)))
#endif

namespace ReVoice
{
  void *rvAlloc(size_t size, const char *file, const char *func, int line);
  void rvFree(void *ptr);
  void rvCheckAllocated();

  template<typename T>static inline T *rvAlloc(size_t n, const char *file, const char *func, int line)
  { return reinterpret_cast<T*>(rvAlloc(n * sizeof(T), file, func, line)); }

  template<typename T>static inline void rvFree(T *&ptr)
  {
    rvFree(reinterpret_cast<void*>(ptr));
    ptr = nullptr;
  }

  template<typename T>static inline T clip(T min, T v, T max)
  { return std::min(std::max(min, v), max); }

  template<typename T>static inline T roundUpToPowerOf2(T v)
  {
    return std::pow(2.0, std::ceil(std::log2(v)));
  }

  template<>static inline unsigned int roundUpToPowerOf2(unsigned int v)
  {
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++;
    return v;
  }

  template<>static inline unsigned long long roundUpToPowerOf2(unsigned long long v)
  {
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v |= v >> 32;
    v++;
    return v;
  }

  template<>static inline int roundUpToPowerOf2(int v)
  {
    if(v >= 0)
      return static_cast<int>(roundUpToPowerOf2(static_cast<unsigned int>(v)));
    else
      return static_cast<int>(roundUpToPowerOf2(static_cast<double>(v)));
  }

  template<>static inline long long roundUpToPowerOf2(long long v)
  {
    rvAssert(v >= 0);
    return static_cast<long>(roundUpToPowerOf2(static_cast<unsigned long long>(v)));
  }

  template<typename T>static inline T lerp(T x0, T x1, double ratio)
  { return x0 + (x1 - x0) * ratio; }

  template<typename T>static inline T *arange(int begin, int end, T *out = nullptr)
  {
    rvAssert(end > begin, "Bad range.");
    if(!out)
      out = RVALLOC(T, end - begin);
    for(int i = begin; i < end; ++i)
      out[i - begin] = static_cast<T>(i);
    return out;
  }

  template<typename T>static inline T *linspace(T begin, T end, int n, T *out = nullptr)
  {
    rvAssert(n > 0, "n must be greater than 0.");
    if(!out)
      out = RVALLOC(T, n);
    double n_d = static_cast<double>(n - 1);
    double range = end - begin;
    for(int i = 0; i < n; ++i)
      out[i] = static_cast<T>(begin + range * (static_cast<double>(i) / n_d));
    return out;
  }

  template<typename T>static inline int *argsort(const T *array, int n, int *out = nullptr)
  {
    rvAssert(n > 0, "n must be greater than 0.");
    if(!out)
      out = RVALLOC(T, n);
    arange(0, n, out);
    std::sort(out, out + n, [array](T x, T y) {return array[x] < array[y]; });
    return out;
  }

  template<typename T>static inline T max(const T *array, int n)
  {
    return *std::max_element(array, array + n);
  }

  template<typename T>static inline T min(const T *array, int n)
  {
    return *std::min_element(array, array + n);
  }

  template<typename T>static inline int argmax(const T *array, int n)
  {
    return static_cast<int>(std::max_element(array, array + n) - array);
  }

  template<typename T>static inline int argmin(const T *array, int n)
  {
    return static_cast<int>(std::min_element(array, array + n) - array);
  }

  template<typename T>static inline T sum(const T *array, int n, T init = 0)
  {
    return std::accumulate(array, array + n, init);
  }

  template<typename T>static inline T squareSum(const T *array, int n, T init = 0)
  {
    for(int i = 0; i < n; ++i)
      init += array[i] * array[i];
    return init;
  }

  template<typename T>static inline T mean(const T *array, int n, T init = 0)
  { return std::accumulate(array, array + n, init) / static_cast<T>(n); }

  template<typename T>static inline T *getitem(const T *array, int *idx, int nIdx, T *out = nullptr)
  {
    rvAssert(nIdx > 0, "nIdx must be greater than 0.");
    if(!out)
      out = RVALLOC(nIdx, end - begin);
    for(int i = 0; i < nIdx; ++i)
      out[i] = array[idx[i]];
    return out;
  }

  static inline RvReal sinc(RvReal x)
  {
    if(x == 0.0)
      return 1.0;
    auto px = 3.1415926535897932384626433832795 * x;
    return std::sin(px) / px;
  }

  template<typename T>static inline void printArray(const T *array, int n)
  {
    rvAssert(n >= 0, "n cannot be less than 0.");
    std::cout << "[";
    for(int i = 0; i < n; ++i)
    {
      if(i != n - 1)
        std::cout << array[i] << ", ";
      else
        std::cout << array[i];
    }
    std::cout << "]" << std::endl;
  }

  template<typename T>static inline bool isInPlace(const T *a, int nA, const T *b)
  { return b >= a && b < a + nA; }

  enum MessageCategory
  {
    DebugMessage = 0,
    WarningMessage,
    CriticalMessage,
    FatalMessage
  };

  void debug(const char *msg, ...);
  void warning(const char *msg, ...);
  void critical(const char *msg, ...);
  void fatal(const char *msg, ...);

  void vPrintMessage(MessageCategory category, const char *msg, va_list args);
} // namespace ReVoice