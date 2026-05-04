#pragma once

#include <cstdint>
#include <string>

namespace bamsi {

struct VersionInfo {
    std::string version;
    std::uint32_t format_version;
};

}  // namespace bamsi
