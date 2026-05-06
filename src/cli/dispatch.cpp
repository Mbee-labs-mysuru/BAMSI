#include "bamsi/cli/dispatch.hpp"

#include <iostream>
#include <string>
#include <unordered_set>

namespace bamsi {

Status DispatchKnownCommand(const std::string& cmd) {
    static const std::unordered_set<std::string> kKnownCommands = {
        "build",
        "count",
        "exists",
        "locate",
        "region-count",
        "region-exists",
        "reconstruct",
        "info",
        "verify",
    };

    if (kKnownCommands.contains(cmd)) {
        std::cerr << "error: subcommand not implemented yet: " << cmd << "\n";
        return Status(StatusCode::kInvalidArgument,
                      "subcommand not implemented yet: " + cmd);
    }

    std::cerr << "error: unknown subcommand: " << cmd << "\n";
    std::cerr << "usage: bamsi <subcommand>\n";
    std::cerr << "subcommands: version, build, count, exists, locate, region-count, region-exists, reconstruct, info, verify\n";
    return Status(StatusCode::kInvalidArgument, "unknown subcommand: " + cmd);
}

}  // namespace bamsi
