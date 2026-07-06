/**
 * fuzz_format_parser.cpp — libFuzzer harness for .bsi format parser
 *
 * Contract §8.3.1: "Run AFL or libFuzzer against the format parser for
 * 7 days continuous on a dedicated machine."
 *
 * This harness feeds random byte sequences to ReadBsi() to find crashes,
 * buffer overflows, and assertion failures in the binary parser.
 *
 * Build:  cmake -DBUILD_FUZZ=ON .. && make fuzz_format_parser
 * Run:    ./fuzz_format_parser corpus/ -max_len=65536 -timeout=10
 *
 * For the 7-day campaign: see scripts/run_fuzzer.sh
 */

#include "format/format.hpp"
#include "bamsix/types.hpp"

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <filesystem>

namespace fs = std::filesystem;

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
    // Write the fuzz input to a temporary file since ReadBsi expects a path
    static int counter = 0;
    std::string tmp_path = "/tmp/bamsix_fuzz_" + std::to_string(getpid()) + "_" +
                           std::to_string(counter++) + ".bsi";

    {
        std::ofstream ofs(tmp_path, std::ios::binary);
        if (!ofs) return 0;
        ofs.write(reinterpret_cast<const char*>(data), size);
    }

    try {
        auto idx = bamsix::ReadBsi(tmp_path);
        // If we get here, the file parsed successfully — no crash
        (void)idx;
    } catch (const bamsix::Error&) {
        // Expected: malformed input should throw Error, not crash
    } catch (const std::exception&) {
        // Also acceptable: standard exceptions
    } catch (...) {
        // Any other exception is also acceptable (no crash = success)
    }

    // Clean up
    fs::remove(tmp_path);
    return 0;
}
