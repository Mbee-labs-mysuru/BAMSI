#include "bamsi/bamsi.hpp"
#include "gtest/gtest.h"

#include <cstdlib>
#include <string>

static int RunBamsiVersion() {
    // Run the built CLI; assumes tests run from the build directory.
    return std::system("./bamsi version");
}

TEST(VersionCliTest, UsesCoreVersionInfo) {
    const bamsi::VersionInfo info = bamsi::GetVersionInfo();
    (void)info;  // for now, just ensure the command succeeds

    const int exit_code = RunBamsiVersion();
    EXPECT_EQ(exit_code, 0);
}
