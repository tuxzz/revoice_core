#include "util_p.hpp"
#include <cstdlib>
#include <cstdarg>

using namespace ReVoice;

namespace ReVoice
{
  void debug(const char * msg, ...)
  {
    va_list args;
    va_start(args, msg);
    vPrintMessage(DebugMessage, msg, args);
    va_end(args);
  }

  void warning(const char * msg, ...)
  {
    va_list args;
    va_start(args, msg);
    vPrintMessage(WarningMessage, msg, args);
    va_end(args);
  }

  void critical(const char * msg, ...)
  {
    va_list args;
    va_start(args, msg);
    vPrintMessage(CriticalMessage, msg, args);
    va_end(args);
  }

  void fatal(const char *msg, ...)
  {
    va_list args;
    va_start(args, msg);
    vPrintMessage(FatalMessage, msg, args);
    va_end(args);
  }

  void vPrintMessage(MessageCategory category, const char *msg, va_list args)
  {
    char buffer[16384] = {'\x00'};
    vsnprintf(buffer, sizeof(buffer), msg, args);
    fprintf(stderr, "%s\n", buffer);

    if(category == FatalMessage)
      std::abort();
  }
} // namespace ReVoice