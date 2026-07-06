#pragma once

#include <string>
#include "bamsix/version.hpp"  // brings in bamsix::VersionInfo

namespace bamsix {

// Do NOT redeclare struct VersionInfo here.
// Use the one from include/bamsix/version.hpp.

VersionInfo GetVersionInfo();      // implementation lives in src/cabi/version.cpp
std::string version_string();      // convenience wrapper

}  // namespace bamsix
