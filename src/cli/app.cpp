#include "bamsi/app.hpp"
#include "bamsi/cli/dispatch.hpp"
#include "bamsi/cli/help.hpp"
#include "bamsi/cli/version.hpp"

#include <string>

namespace bamsi {

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

}  // namespace bamsi
