#include "Profile.h"

std::map<const char *, Profiler::ProfileEntry> Profiler::profileMap;
Profiler *Profiler::currentProfiler = NULL;
const char *Profiler::rootName = NULL;
