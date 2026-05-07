#include "bamsi/bamsi.hpp"
#include "bamsi/config.hpp"
#include "cabi/version.hpp"

namespace bamsi {

VersionInfo GetVersionInfo() {
    return VersionInfo{
        .version        = std::string(BAMSI_VERSION),
        .format_version = BAMSI_FORMAT_VERSION,
    };
}

std::string version_string() {
    // Delegate to GetVersionInfo so there is only one source of truth.
    return GetVersionInfo().version;
}

}  // namespace bamsi
