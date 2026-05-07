#pragma once

#include "bamsi/types.hpp"
#include "../ingest/ingest.hpp"
#include <vector>

namespace bamsi {

/// Result of read ordering (Stage 2 of build pipeline).
struct OrderingResult {
    std::vector<OrderedRead>    reads;
    std::array<uint8_t, 32>     ordering_hash;
};

/// Sort reads by (chrom_id, pos, source_file_id, bam_offset), assign read_id,
/// compute ordering_hash per §0.10.
OrderingResult OrderReads(const IngestResult& ingest);

}  // namespace bamsi
