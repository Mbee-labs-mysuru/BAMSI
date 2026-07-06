#pragma once

#include "bamsix/types.hpp"

#include <string>

namespace bamsix {

struct VersionInfo {
    std::string version;
    std::uint32_t format_version;
};

VersionInfo GetVersionInfo();

}  // namespace bamsix
