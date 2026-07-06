#include "bamsix/bamsix.hpp"
#include "bamsix/app.hpp"
#include <iostream>

namespace bamsix {

Status RunVersion() {
    const auto info = GetVersionInfo();
    std::cout << "bamsix " << info.version
              << " format-version " << info.format_version << "\n";
    return Status::Ok();
}

}  // namespace bamsix
