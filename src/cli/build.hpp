#pragma once

#include "bamsi/types.hpp"
#include <string>
#include <vector>

namespace bamsi {

/// Build a .bsi index from BAM input files.
/// This is the top-level orchestrator for the full build pipeline
/// (Stages 1–10 of Architecture §4).
void BuildIndex(const std::vector<std::string>& bam_paths,
                const std::string& output_path,
                const BuildConfig& config);

}  // namespace bamsi
