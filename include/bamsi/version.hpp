#pragma once

#include "bamsi/types.hpp"

#include <string>

namespace bamsi {

struct VersionInfo {
    std::string version;
    std::uint32_t format_version;
};

VersionInfo GetVersionInfo();

}  // namespace bamsi
