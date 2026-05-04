#include "bamsi/bamsi.hpp"
#include "bamsi/app.hpp"
#include <iostream>

namespace bamsi {

Status RunVersion() {
    const auto info = GetVersionInfo();
    std::cout << "bamsi " << info.version
              << " format-version " << info.format_version << "\n";
    return Status::Ok();
}

}  // namespace bamsi
