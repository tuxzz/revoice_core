#include "util_p.hpp"
#include <cstdlib>
#include <unordered_map>
#include <mutex>

#ifdef _MSC_VER
#pragma warning(disable: 4127)
#endif

using namespace ReVoice;

constexpr static size_t alignment = 64;

struct MemoryBlockAfter
{
  size_t magic;
};

struct MemoryBlockBefore
{
  void *headPtr;
  const char *file, *func;
  int line;
  MemoryBlockAfter *blockAfter;
  size_t size;
  size_t magic;
};

#ifndef DISABLE_ALLOCGUARD
static std::mutex g_memoryBlockDictLock;
static std::unordered_map<void*, MemoryBlockBefore*> g_memoryBlockDict;
#endif // DISABLE_ALLOCGUARD

constexpr static size_t magicNumber()
{
  static_assert(sizeof(size_t) == 8 || sizeof(size_t) == 4, "Unsupported size_t size.");
  if(sizeof(size_t) == 8)
    return 0x1234567852391314;
  else
    return 0x52391314;
}

static void verifyMemoryBlock(MemoryBlockBefore *before, void *ptr)
{
  if(before->magic != magicNumber() || reinterpret_cast<char*>(ptr) + before->size != reinterpret_cast<char*>(before->blockAfter) || before->blockAfter->magic != magicNumber())
  {
    fprintf(stderr, "\nRVD: Memory block 0x%p was broken!\n", ptr);
    if(before->file && before->func)
      fprintf(stderr, "     Memory block is allocated at %s@%s:%d\n", before->file, before->func, before->line);
    else if(before->file)
      fprintf(stderr, "     Memory block is allocated at %s\n", before->file);
    else
      fprintf(stderr, "     Memory block allocation position is not clear.");
    std::abort();
  }
}

namespace ReVoice
{
  void *rvAlloc(size_t size, const char *file, const char *func, int line)
  {
    rvAssert(size > 0, "size must be greater than 0");
    size_t additionalSize = sizeof(MemoryBlockBefore) + sizeof(MemoryBlockAfter);
    size_t paddedSize = size + alignment - 1;
    size_t totalSize = paddedSize + additionalSize;
    auto headPtr = reinterpret_cast<char*>(malloc(totalSize));
    auto additionalPtr = headPtr + sizeof(MemoryBlockBefore);

    auto alignedAdditionalPtr = additionalPtr + alignment - (reinterpret_cast<size_t>(additionalPtr) % alignment);
    auto before = reinterpret_cast<MemoryBlockBefore*>(alignedAdditionalPtr - sizeof(MemoryBlockBefore));
    auto after = reinterpret_cast<MemoryBlockAfter*>(alignedAdditionalPtr + size);

    before->headPtr = headPtr;
    before->file = file;
    before->func = func;
    before->line = line;
    before->blockAfter = after;
    before->size = size;
    before->magic = magicNumber();
    after->magic = magicNumber();

#ifndef DISABLE_ALLOCGUARD
    std::unique_lock<std::mutex> locker(g_memoryBlockDictLock);
    g_memoryBlockDict.insert(std::make_pair(alignedAdditionalPtr, before));
    locker.unlock();
#endif // DISABLE_ALLOCGUARD

    return reinterpret_cast<void*>(alignedAdditionalPtr);
  }

  void rvFree(void *ptr)
  {
    rvAssert(ptr, "cannot free a null pointer");
    auto alignedAdditionalPtr = reinterpret_cast<char*>(ptr);
    auto before = reinterpret_cast<MemoryBlockBefore*>(alignedAdditionalPtr - sizeof(MemoryBlockBefore));
    verifyMemoryBlock(before, ptr);

#ifndef DISABLE_ALLOCGUARD
    std::unique_lock<std::mutex> locker(g_memoryBlockDictLock);
    g_memoryBlockDict.erase(g_memoryBlockDict.find(ptr));
    locker.unlock();
#endif // DISABLE_ALLOCGUARD

    free(before->headPtr);
  }

#ifndef DISABLE_ALLOCGUARD
  void rvCheckAllocated()
  {
    std::unique_lock<std::mutex> locker(g_memoryBlockDictLock);
    for(auto &pair : g_memoryBlockDict)
    {
      void *ptr = pair.first;
      MemoryBlockBefore *before = pair.second;
      verifyMemoryBlock(before, ptr);
      if(before->file && before->func)
        critical("Memory block %p(allocated at %s@%s:%d) with size %lu", ptr, before->file, before->func, before->line, before->size);
      else if(before->file)
        critical("Memory block %p(allocated at %s) with size %lu", ptr, before->file, before->size);
      else
        critical("Memory block %p(allocated at unknown position) with size %lu", ptr);
    }
    locker.unlock();
  }
#else // DISABLE_ALLOCGUARD
  void rvCheckAllocated()
  {}
#endif // DISABLE_ALLOCGUARD
} // namespace ReVoice