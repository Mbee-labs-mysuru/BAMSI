#include "bamsix/bamsix.hpp"
#include "gtest/gtest.h"
#include <regex>

namespace {

TEST(VersionInfoTest, HasSemanticVersionAndFormat) {
    const bamsix::VersionInfo info = bamsix::GetVersionInfo();

    // version is non-empty and looks like MAJOR.MINOR.PATCH
    EXPECT_FALSE(info.version.empty());

    const std::regex semantic_version_pattern{R"(^\d+\.\d+\.\d+(-.+)?$)"};
    EXPECT_TRUE(std::regex_match(info.version, semantic_version_pattern))
        << "version string was: " << info.version;

    // format_version is a positive integer (tighten later if needed)
    EXPECT_GT(info.format_version, 0);
}

}  // namespace
