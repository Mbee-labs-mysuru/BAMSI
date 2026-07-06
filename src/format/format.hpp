#pragma once

#include "bamsix/types.hpp"
#include "../fmindex/fmindex.hpp"
#include "../bitvectors/bitvectors.hpp"
#include "../sarange/sarange.hpp"
#include <map>
#include <string>
#include <vector>

namespace bamsix {

/// All data loaded from a .bsi file, sufficient for query operations.
struct LoadedIndex {
    BsiHeader                           header;
    FMIndexEngine                       fm;
    BitVectors                          bv;
    WindowTable                         windows;
    std::vector<std::string>            chrom_names;
    std::map<std::string, uint32_t>     chrom_to_id;

    // Read metadata (for Locate/RegionalCount — chrom_id, pos, cigar, seq_len)
    std::vector<OrderedRead>            reads;

    // Stream payloads (for reconstruction)
    std::vector<uint8_t>                seq_payload;
    std::vector<uint8_t>                qual_payload;
    std::vector<uint8_t>                meta_payload;
    std::vector<uint8_t>                map_payload;

    // Per-read stream directories (C5 fix: loaded from .bsi for reconstruction)
    // Contract §2.3 F4: qual may be per-read or block-level
    StreamDirectory                     qual_directory;
    uint32_t                            qual_block_size = 0;  // 0 = per-read
    StreamDirectoryPerRead              meta_directory;
    StreamDirectoryPerRead              map_directory;

    // ENHANCED tier: SARange wavelet tree (Architecture §5.3)
    SARange                             sarange;

    // Checksum validation results
    bool                                global_checksum_valid = false;
};

/// Read and parse a .bsi file. Throws Error on parse failure.
LoadedIndex ReadBsi(const std::string& path);

/// Verify a .bsi file's integrity (checksums only, no full load).
/// Returns true if all checksums pass.
bool VerifyBsi(const std::string& path);

}  // namespace bamsix
