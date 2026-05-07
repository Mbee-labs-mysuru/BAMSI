#pragma once

#include "bamsi/app.hpp"
#include <string>

namespace bamsi {

Status DispatchKnownCommand(const std::string& cmd);
Status DispatchCommand(const std::string& cmd, int argc, char** argv);

}  // namespace bamsi
