#include "bamsix/app.hpp"
#include "bamsix/cli/dispatch.hpp"
#include "bamsix/cli/help.hpp"
#include "bamsix/cli/version.hpp"

#include <string>

namespace bamsix {

Status RunApp(int argc, char** argv) {
    if (argc <= 1) {
        return RunHelp();
    }

    const std::string cmd = argv[1];

    if (cmd == "--help" || cmd == "-h" || cmd == "help") {
        return RunHelp();
    }

    if (cmd == "version") {
        return RunVersion();
    }

    return DispatchCommand(cmd, argc, argv);
}

}  // namespace bamsix
