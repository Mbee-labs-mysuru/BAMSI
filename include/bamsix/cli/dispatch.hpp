#pragma once

#include "bamsix/app.hpp"
#include <string>

namespace bamsix {

Status DispatchKnownCommand(const std::string& cmd);
Status DispatchCommand(const std::string& cmd, int argc, char** argv);

}  // namespace bamsix
