#include "bamsix/app.hpp"
#include <iostream>

namespace bamsix {

Status RunHelp() {
    std::cout << "usage: bamsix <subcommand>\n";
    std::cout << "subcommands: version, build, count, exists, locate, region-count, region-exists, reconstruct, info, verify\n";
    return Status::Ok();
}

}  // namespace bamsix
