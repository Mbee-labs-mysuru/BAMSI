#include "mapping.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>

namespace bamsi {

// ─── Separator detection via bitvector only (I14) ────────────────────────────

bool IsSeparatorPosition(uint64_t pos, const SuccinctBitvector& B_read) {
    // Architecture §5.1: A separator '#' is at position readStarts[i]-1 for i≥1.
    // Bitvector-only check (I14 — no raw S access):
    //   pos is a separator iff:
    //   1. pos > 0 (position 0 is always the start of read 0)
    //   2. B_read[pos] == 0 (pos is NOT a read start)
    //   3. pos + 1 < |B_read| AND B_read[pos+1] == 1
    //      (the next position IS a read start, so pos is the '#' before it)
    if (pos == 0) return false;
    if (B_read.Access(pos)) return false;  // pos is a readStart, not a separator
    if (pos + 1 < B_read.Length() && B_read.Access(pos + 1)) return true;
    return false;
}

// ─── CIGAR mapping (Architecture §5.2) ──────────────────────────────────────

uint64_t CigarRefPos(const CigarRecord& cigar, uint64_t p_anchor,
                     uint64_t offset, CigarDirection direction) {
    // Traverse CIGAR ops maintaining ref_pos and read_pos.
    uint64_t ref_pos  = p_anchor;
    uint64_t read_pos = 0;

    // Track aligned bases for fallback
    uint64_t last_aligned_ref = 0;
    bool     has_last_aligned = false;
    uint64_t first_aligned_ref_after = 0;
    bool     found_aligned_after = false;

    // First pass: check if there are ANY aligned bases
    bool has_any_aligned = false;
    for (const auto& op : cigar) {
        if (op.op == 0 || op.op == 7 || op.op == 8) {  // M, =, X
            has_any_aligned = true;
            break;
        }
    }

    if (!has_any_aligned) {
        // No aligned base: return p_anchor per §5.2 fallback
        return p_anchor;
    }

    // Main traversal
    ref_pos  = p_anchor;
    read_pos = 0;
    last_aligned_ref = p_anchor;
    has_last_aligned = false;

    for (size_t ci = 0; ci < cigar.size(); ++ci) {
        const auto& cop = cigar[ci];

        switch (cop.op) {
            case 0:  // M
            case 7:  // =
            case 8:  // X
            {
                // Consumes both read and reference
                if (offset >= read_pos && offset < read_pos + cop.len) {
                    return ref_pos + (offset - read_pos);
                }
                last_aligned_ref = ref_pos + cop.len - 1;
                has_last_aligned = true;
                ref_pos  += cop.len;
                read_pos += cop.len;
                break;
            }
            case 1:  // I (insertion to reference)
            {
                if (offset >= read_pos && offset < read_pos + cop.len) {
                    if (direction == CigarDirection::LEFT) {
                        // Return last aligned ref base before this I
                        return has_last_aligned ? last_aligned_ref : p_anchor;
                    } else {
                        // Return first aligned ref base after this I
                        // Look ahead for next aligned op
                        uint64_t look_ref = ref_pos;
                        uint64_t look_read = read_pos + cop.len;
                        for (size_t k = ci + 1; k < cigar.size(); ++k) {
                            if (cigar[k].op == 0 || cigar[k].op == 7 || cigar[k].op == 8) {
                                return look_ref;
                            }
                            if (cigar[k].op == 2 || cigar[k].op == 3) {
                                look_ref += cigar[k].len;
                            } else if (cigar[k].op == 1) {
                                look_read += cigar[k].len;
                            } else if (cigar[k].op == 4) {
                                look_read += cigar[k].len;
                            }
                            // H, P: no-op
                        }
                        // No aligned base after: fall back to last before
                        return has_last_aligned ? last_aligned_ref : p_anchor;
                    }
                }
                read_pos += cop.len;
                break;
            }
            case 2:  // D
            case 3:  // N
            {
                ref_pos += cop.len;
                break;
            }
            case 4:  // S (soft clip)
            {
                if (offset >= read_pos && offset < read_pos + cop.len) {
                    // Map to nearest aligned reference base
                    // Find nearest aligned base before and after
                    uint64_t best_ref = p_anchor;
                    uint64_t best_dist = UINT64_MAX;

                    // Last aligned before
                    if (has_last_aligned) {
                        uint64_t dist = offset - read_pos;  // approximate
                        if (dist < best_dist || (dist == best_dist && last_aligned_ref < best_ref)) {
                            best_dist = dist;
                            best_ref = last_aligned_ref;
                        }
                    }

                    // First aligned after: look ahead
                    uint64_t look_ref = ref_pos;
                    uint64_t look_read = read_pos + cop.len;
                    for (size_t k = ci + 1; k < cigar.size(); ++k) {
                        if (cigar[k].op == 0 || cigar[k].op == 7 || cigar[k].op == 8) {
                            uint64_t dist_after = look_read - offset;
                            if (dist_after < best_dist ||
                                (dist_after == best_dist && look_ref < best_ref)) {
                                best_ref = look_ref;
                            }
                            break;
                        }
                        if (cigar[k].op == 2 || cigar[k].op == 3) {
                            look_ref += cigar[k].len;
                        } else if (cigar[k].op == 1 || cigar[k].op == 4) {
                            look_read += cigar[k].len;
                        }
                    }

                    return best_ref;
                }
                read_pos += cop.len;
                break;
            }
            case 5:  // H (hard clip)
            case 6:  // P (padding)
            {
                // No-op: neither pointer changes
                break;
            }
            default:
                break;
        }
    }

    // Should not reach here for valid inputs
    return p_anchor;
}

// ─── MapOccurrence (Architecture §5.1) ──────────────────────────────────────

MappingResult MapOccurrence(uint64_t pos, uint64_t pattern_len,
                            QueryStrand strand,
                            const SuccinctBitvector& B_read,
                            const std::vector<OrderedRead>& reads) {
    // Step 1: Identify read (closed-interval rank convention)
    uint64_t read_id = B_read.Rank1(pos) - 1;

    // Step 2: Read start and match offset
    uint64_t read_start   = B_read.Select1(read_id + 1);
    uint64_t offset_start = pos - read_start;
    uint64_t offset_end   = offset_start + pattern_len - 1;

    // Step 3: Retrieve coordinates from S_map (chrom_id, pos)
    const auto& read = reads[read_id];
    uint32_t chrom_id = read.chrom_id;
    uint64_t p_i = read.pos;  // 1-based SAM POS

    // Step 4: Retrieve CIGAR from S_meta
    const auto& cigar = read.cigar;

    // Step 5: Apply CIGAR mapping
    uint64_t p_min = CigarRefPos(cigar, p_i, offset_start, CigarDirection::LEFT);
    uint64_t p_max = CigarRefPos(cigar, p_i, offset_end, CigarDirection::RIGHT);

    // Ensure p_min <= p_max
    if (p_min > p_max) std::swap(p_min, p_max);

    return MappingResult{
        .chrom_id     = chrom_id,
        .p_min        = p_min,
        .p_max        = p_max,
        .read_id      = read_id,
        .query_strand = strand,
    };
}

}  // namespace bamsi
