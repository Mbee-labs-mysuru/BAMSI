#include "bamsi/bamsi.hpp"
#include "bamsi/config.hpp"

namespace bamsi {

VersionInfo GetVersionInfo() {
    return VersionInfo{
        .version = std::string(BAMSI_VERSION),
        .format_version = BAMSI_FORMAT_VERSION,
    };
}

}  // namespace bamsi
