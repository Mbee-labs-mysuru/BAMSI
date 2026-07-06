#pragma once

#include "bamsix/types.hpp"
#include <string>
#include <vector>
#include <map>

namespace bamsix {

/// Result of BAM ingestion (Stage 1 of build pipeline).
struct IngestResult {
    std::vector<RawRead>          reads;
    std::map<std::string, uint32_t> chrom_to_id;   // chrom name → chrom_id
    std::vector<std::string>      chrom_names;     // chrom_id → chrom name (sorted)
    uint32_t                      source_file_count;
    std::array<uint8_t, 32>       source_manifest_hash;
};

/// Ingest one or more BAM files.
/// Applies the inclusion rule (§0.1): mapped, primary, non-supplementary, POS >= 1.
/// Normalises bases to codes 0..4, converts BAM 0-based POS to 1-based SAM POS.
/// @param bam_paths  Ordered list of BAM file paths (order defines source_file_id).
/// @return IngestResult containing all reads in ℛ and chromosome mapping.
IngestResult IngestBams(const std::vector<std::string>& bam_paths);

}  // namespace bamsix
