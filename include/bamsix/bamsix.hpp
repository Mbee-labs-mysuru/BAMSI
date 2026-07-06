#pragma once

#include "bamsix/format.hpp"
#include "bamsix/result.hpp"
#include "bamsix/status.hpp"
#include "bamsix/types.hpp"
#include "bamsix/version.hpp"

namespace bamsix {
	VersionInfo GetVersionInfo();
	Status RunCli(int argc, char** argv);
}  // namespace bamsix
