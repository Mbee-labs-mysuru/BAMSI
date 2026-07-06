#include "bamsix/bamsix.hpp"
#include "bamsix/config.hpp"
#include "cabi/version.hpp"

namespace bamsix {

VersionInfo GetVersionInfo() {
    return VersionInfo{
        .version        = std::string(BAMSIX_VERSION),
        .format_version = BAMSIX_FORMAT_VERSION,
    };
}

std::string version_string() {
    // Delegate to GetVersionInfo so there is only one source of truth.
    return GetVersionInfo().version;
}

}  // namespace bamsix
