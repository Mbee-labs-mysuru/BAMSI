#pragma once

#include "bamsi/format.hpp"
#include "bamsi/result.hpp"
#include "bamsi/status.hpp"
#include "bamsi/types.hpp"

namespace bamsi {

VersionInfo GetVersionInfo();
Status RunCli(int argc, char** argv);

}  // namespace bamsi
