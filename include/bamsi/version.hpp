#pragma once

#include "bamsi/types.hpp"

namespace bamsi {

struct VersionInfo {
    std::string version;
    uint32_t format_version;
};

VersionInfo GetVersionInfo();

}  // namespace bamsi
