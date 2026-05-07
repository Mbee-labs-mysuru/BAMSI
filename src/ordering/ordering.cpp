#include "ordering.hpp"

#include <algorithm>
#include <openssl/sha.h>

namespace bamsi {

OrderingResult OrderReads(const IngestResult& ingest) {
    OrderingResult result;

    // Copy reads and assign chrom_id
    result.reads.reserve(ingest.reads.size());
    for (const auto& raw : ingest.reads) {
        OrderedRead ordered;
        // Copy all RawRead fields
        ordered.seq            = raw.seq;
        ordered.chrom          = raw.chrom;
        ordered.flag           = raw.flag;
        ordered.pos            = raw.pos;
        ordered.cigar          = raw.cigar;
        ordered.qual           = raw.qual;
        ordered.source_file_id = raw.source_file_id;
        ordered.bam_offset     = raw.bam_offset;
        // Assign chrom_id from chrom_to_id map
        auto it = ingest.chrom_to_id.find(raw.chrom);
        if (it != ingest.chrom_to_id.end()) {
            ordered.chrom_id = it->second;
        } else {
            ordered.chrom_id = 0;  // Should not happen
        }
        ordered.read_id = 0; // Will be assigned after sorting
        result.reads.push_back(std::move(ordered));
    }

    // Sort by (chrom_id, pos, source_file_id, bam_offset) per §2.4
    std::stable_sort(result.reads.begin(), result.reads.end(),
        [](const OrderedRead& a, const OrderedRead& b) {
            if (a.chrom_id != b.chrom_id) return a.chrom_id < b.chrom_id;
            if (a.pos != b.pos) return a.pos < b.pos;
            if (a.source_file_id != b.source_file_id)
                return a.source_file_id < b.source_file_id;
            return a.bam_offset < b.bam_offset;
        });

    // Assign read_id = 0-based rank in sorted order
    for (uint64_t i = 0; i < result.reads.size(); ++i) {
        result.reads[i].read_id = i;
    }

    // Compute ordering_hash per Contract §0.10 / Architecture §2.4:
    // SHA-256(concat for i=0..N-1 in read_id order:
    //   uint32_le(chrom_id) || uint64_le(pos) ||
    //   uint32_le(source_file_id) || uint64_le(bam_offset))
    // NOTE: read_id is NOT part of the hash — it is implicit in the ordering.
    {
        SHA256_CTX ctx;
        SHA256_Init(&ctx);
        for (const auto& r : result.reads) {
            uint32_t cid  = r.chrom_id;
            uint64_t pos  = r.pos;
            uint32_t sfid = r.source_file_id;
            uint64_t boff = r.bam_offset;
            SHA256_Update(&ctx, &cid,  sizeof(cid));
            SHA256_Update(&ctx, &pos,  sizeof(pos));
            SHA256_Update(&ctx, &sfid, sizeof(sfid));
            SHA256_Update(&ctx, &boff, sizeof(boff));
        }
        SHA256_Final(result.ordering_hash.data(), &ctx);
    }

    return result;
}

}  // namespace bamsi
