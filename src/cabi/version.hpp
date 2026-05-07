#pragma once

#include <string>
#include "bamsi/version.hpp"  // brings in bamsi::VersionInfo

namespace bamsi {

// Do NOT redeclare struct VersionInfo here.
// Use the one from include/bamsi/version.hpp.

VersionInfo GetVersionInfo();      // implementation lives in src/cabi/version.cpp
std::string version_string();      // convenience wrapper

}  // namespace bamsi
