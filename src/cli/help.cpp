#include "bamsi/app.hpp"
#include <iostream>

namespace bamsi {

Status RunHelp() {
    std::cout << "usage: bamsi <subcommand>\n";
    std::cout << "subcommands: version, build, count, exists, locate, region-count, region-exists, reconstruct, info, verify\n";
    return Status::Ok();
}

}  // namespace bamsi
