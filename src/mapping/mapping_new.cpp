#include "mapping.hpp"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>

namespace bamsix {

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

// ─── BuildAbsolutePosTable ────────────────────────────────────────────────────
// Called ONCE at query load time. Scans all N map blocks in order and records
// one entry per chromosome boundary (is_delta==false). Cost: O(N) one-time.
// After this, each MapOccurrenceLazy call costs O(log C + K) where C is the
// number of chromosomes and K is reads between two consecutive boundaries.

std::vector<AbsPosEntry> BuildAbsolutePosTable(
    const std::vector<uint8_t>& map_payload,
    const StreamDirectoryPerRead& dir_map,
    uint8_t map_codec_id) {

    std::vector<AbsPosEntry> table;
    uint64_t N = dir_map.size();
    if (N == 0) return table;

    // Reserve conservatively (at most one entry per chromosome; human has ~25)
    table.reserve(64);

    for (uint64_t i = 0; i < N; ++i) {
        DecodedMap dm = DecodeMapRead(map_payload, dir_map[i], map_codec_id);
        if (!dm.is_delta) {
            // This read stores an absolute position — record it
            table.push_back({i, dm.pos});
        }
    }

    // Table is already sorted by read_id (we iterate in order)
    return table;
}

// ─── BuildFullAbsPosTable ─────────────────────────────────────────────────────
// Single O(N) sequential pass: decodes every map block, accumulates deltas,
// stores the resulting absolute position for each read.
// Memory: N × 8 bytes (16 MB for 2M reads — negligible).
// After this, MapOccurrenceLazy does positions[read_id] — O(1), no loops.

FullAbsPosTable BuildFullAbsPosTable(
    const std::vector<uint8_t>& map_payload,
    const StreamDirectoryPerRead& dir_map,
    uint8_t map_codec_id) {

    FullAbsPosTable table;
    uint64_t N = dir_map.size();
    if (N == 0) return table;

    table.positions.resize(N, 0);
    uint64_t acc = 0;  // running absolute position

    for (uint64_t i = 0; i < N; ++i) {
        DecodedMap dm = DecodeMapRead(map_payload, dir_map[i], map_codec_id);
        if (!dm.is_delta) {
            // Absolute entry (chromosome boundary): reset accumulator
            acc = dm.pos;
        } else {
            // Delta entry: accumulate
            int64_t delta;
            std::memcpy(&delta, &dm.pos, sizeof(delta));
            acc = static_cast<uint64_t>(static_cast<int64_t>(acc) + delta);
        }
        table.positions[i] = acc;
    }

    return table;
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

// ─── MapOccurrenceLazy (Contract §2.3 Table B / §3.6.2 Steps 3-4) ──────────
//
// Delta resolution strategy:
//
//   WITH full_abs_table (pre-built by BuildFullAbsPosTable at load time):
//     Direct O(1) array lookup: p_i = full_abs_table->positions[read_id]
//     No loops, no decompression, no delta accumulation.
//
//   WITHOUT full_abs_table (legacy fallback):
//     Performs original O(N) backward+forward scan. Preserved for compatibility.

MappingResult MapOccurrenceLazy(uint64_t pos, uint64_t pattern_len,
                                QueryStrand strand,
                                const SuccinctBitvector& B_read,
                                const std::vector<uint8_t>& meta_payload,
                                const std::vector<uint8_t>& map_payload,
                                const StreamDirectoryPerRead& dir_meta,
                                const StreamDirectoryPerRead& dir_map,
                                uint8_t meta_codec_id,
                                uint8_t map_codec_id,
                                const FullAbsPosTable* full_abs_table) {
    // Step 1: Identify read via bitvector (closed-interval rank convention)
    uint64_t read_id = B_read.Rank1(pos) - 1;

    // Step 2: Read start and match offset
    uint64_t read_start   = B_read.Select1(read_id + 1);
    uint64_t offset_start = pos - read_start;
    uint64_t offset_end   = offset_start + pattern_len - 1;

    // Step 3: Retrieve coordinates from S_map
    DecodedMap map_data = DecodeMapRead(map_payload, dir_map[read_id], map_codec_id);

    uint32_t chrom_id = map_data.chrom_id;
    uint64_t p_i = map_data.pos;

    if (map_data.is_delta) {
        if (full_abs_table && read_id < full_abs_table->positions.size()) {
            // ── O(1) fast path: direct lookup from pre-computed dense table ──
            p_i = full_abs_table->positions[read_id];
        } else {
            // ── Legacy fallback: O(N) backward + forward scan (no table) ──
            uint64_t abs_idx = read_id;
            uint64_t abs_pos = 0;
            while (abs_idx > 0) {
                --abs_idx;
                DecodedMap prev = DecodeMapRead(map_payload, dir_map[abs_idx], map_codec_id);
                if (!prev.is_delta) {
                    abs_pos = prev.pos;
                    break;
                }
            }
            if (abs_idx == 0) {
                DecodedMap first = DecodeMapRead(map_payload, dir_map[0], map_codec_id);
                if (!first.is_delta) {
                    abs_pos = first.pos;
                }
            }
            uint64_t acc = abs_pos;
            for (uint64_t j = abs_idx + 1; j <= read_id; ++j) {
                DecodedMap dm = DecodeMapRead(map_payload, dir_map[j], map_codec_id);
                if (!dm.is_delta) {
                    acc = dm.pos;
                } else {
                    int64_t delta;
                    std::memcpy(&delta, &dm.pos, sizeof(delta));
                    acc = static_cast<uint64_t>(static_cast<int64_t>(acc) + delta);
                }
            }
            p_i = acc;
        }
    }

    // Step 4: Retrieve CIGAR from S_meta (Contract §3.6.2 Step 4)
    //   "cigar ← decode per-read block from S_meta using dir_meta[read_id]"
    DecodedMeta meta_data = DecodeMetaRead(meta_payload, dir_meta[read_id], meta_codec_id);

    // Step 5: Apply CIGAR mapping
    uint64_t p_min = CigarRefPos(meta_data.cigar, p_i, offset_start, CigarDirection::LEFT);
    uint64_t p_max = CigarRefPos(meta_data.cigar, p_i, offset_end, CigarDirection::RIGHT);

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

}  // namespace bamsix

