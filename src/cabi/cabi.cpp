/// BAMSI C ABI Implementation — Contract §10.3
/// All functions are extern "C", thread-safe for concurrent reads.

#include "bamsi/bamsi.h"
#include "bamsi/config.hpp"
#include "format/format.hpp"
#include "query/query.hpp"
#include "mapping/mapping.hpp"

#include <cstring>
#include <memory>
#include <string>
#include <vector>

using namespace bamsi;

/// Opaque index handle — wraps a LoadedIndex.
struct bamsi_index {
    LoadedIndex idx;
};

// ─── Version ───────────────────────────────────────────────────────────

extern "C" const char* bamsi_version(void) {
    return BAMSI_VERSION;
}

extern "C" uint16_t bamsi_format_version(void) {
    return BAMSI_FORMAT_VERSION;
}

// ─── Index lifecycle ───────────────────────────────────────────────────

extern "C" bamsi_status_t bamsi_open(const char* path, bamsi_index_t** out) {
    if (!path || !out) return BAMSI_STATUS_INVALID_ARGUMENT;
    try {
        auto* handle = new bamsi_index();
        handle->idx = ReadBsi(std::string(path));
        *out = handle;
        return BAMSI_STATUS_OK;
    } catch (const Error& e) {
        if (e.code == ErrorCode::CORRUPT_BSI) return BAMSI_STATUS_CORRUPT_BSI;
        return BAMSI_STATUS_INTERNAL_ERROR;
    } catch (...) {
        return BAMSI_STATUS_INTERNAL_ERROR;
    }
}

extern "C" void bamsi_free(bamsi_index_t** idx) {
    if (idx && *idx) {
        delete *idx;
        *idx = nullptr;
    }
}

// ─── Verification ──────────────────────────────────────────────────────

extern "C" bamsi_status_t bamsi_verify(const char* path, int* valid_out) {
    if (!path || !valid_out) return BAMSI_STATUS_INVALID_ARGUMENT;
    try {
        *valid_out = VerifyBsi(std::string(path)) ? 1 : 0;
        return BAMSI_STATUS_OK;
    } catch (...) {
        return BAMSI_STATUS_INTERNAL_ERROR;
    }
}

// ─── Index info ────────────────────────────────────────────────────────

extern "C" bamsi_status_t bamsi_get_n_reads(const bamsi_index_t* idx, uint64_t* out) {
    if (!idx || !out) return BAMSI_STATUS_INVALID_ARGUMENT;
    *out = idx->idx.header.N_reads;
    return BAMSI_STATUS_OK;
}

extern "C" bamsi_status_t bamsi_get_s_length(const bamsi_index_t* idx, uint64_t* out) {
    if (!idx || !out) return BAMSI_STATUS_INVALID_ARGUMENT;
    *out = idx->idx.header.S_length;
    return BAMSI_STATUS_OK;
}

extern "C" bamsi_status_t bamsi_get_n_windows(const bamsi_index_t* idx, uint32_t* out) {
    if (!idx || !out) return BAMSI_STATUS_INVALID_ARGUMENT;
    *out = idx->idx.header.N_windows;
    return BAMSI_STATUS_OK;
}

extern "C" bamsi_status_t bamsi_get_chrom_count(const bamsi_index_t* idx, uint32_t* out) {
    if (!idx || !out) return BAMSI_STATUS_INVALID_ARGUMENT;
    *out = idx->idx.header.chrom_count;
    return BAMSI_STATUS_OK;
}

extern "C" bamsi_status_t bamsi_get_chrom_name(const bamsi_index_t* idx,
                                                uint32_t chrom_idx,
                                                char* buf, size_t buf_len,
                                                size_t* out_len) {
    if (!idx || !out_len) return BAMSI_STATUS_INVALID_ARGUMENT;
    if (chrom_idx >= idx->idx.chrom_names.size()) return BAMSI_STATUS_INVALID_ARGUMENT;
    const auto& name = idx->idx.chrom_names[chrom_idx];
    *out_len = name.size();
    if (buf && buf_len > 0) {
        size_t copy_len = std::min(buf_len - 1, name.size());
        std::memcpy(buf, name.c_str(), copy_len);
        buf[copy_len] = '\0';
    }
    return BAMSI_STATUS_OK;
}

extern "C" bamsi_status_t bamsi_is_lossless(const bamsi_index_t* idx, int* out) {
    if (!idx || !out) return BAMSI_STATUS_INVALID_ARGUMENT;
    *out = idx->idx.header.is_lossless ? 1 : 0;
    return BAMSI_STATUS_OK;
}

// ─── Query: GlobalCount ────────────────────────────────────────────────

extern "C" bamsi_status_t bamsi_global_count(const bamsi_index_t* idx,
                                              const uint8_t* pattern, size_t pat_len,
                                              uint64_t* count) {
    if (!idx || !pattern || pat_len == 0 || !count) return BAMSI_STATUS_INVALID_ARGUMENT;
    try {
        std::vector<uint8_t> pat(pattern, pattern + pat_len);
        StrandMode mode = static_cast<StrandMode>(idx->idx.header.strand_mode);
        *count = GlobalCount(pat, idx->idx.fm, mode);
        return BAMSI_STATUS_OK;
    } catch (...) {
        return BAMSI_STATUS_INTERNAL_ERROR;
    }
}

// ─── Query: GlobalExists ───────────────────────────────────────────────

extern "C" bamsi_status_t bamsi_global_exists(const bamsi_index_t* idx,
                                               const uint8_t* pattern, size_t pat_len,
                                               uint64_t threshold, int* exists) {
    if (!idx || !pattern || pat_len == 0 || !exists) return BAMSI_STATUS_INVALID_ARGUMENT;
    try {
        std::vector<uint8_t> pat(pattern, pattern + pat_len);
        StrandMode mode = static_cast<StrandMode>(idx->idx.header.strand_mode);
        *exists = GlobalExists(pat, idx->idx.fm, mode) ? 1 : 0;
        return BAMSI_STATUS_OK;
    } catch (...) {
        return BAMSI_STATUS_INTERNAL_ERROR;
    }
}

// ─── Query: Locate ─────────────────────────────────────────────────────

extern "C" bamsi_status_t bamsi_locate(const bamsi_index_t* idx,
                                        const uint8_t* pattern, size_t pat_len,
                                        bamsi_locate_result_t* results,
                                        size_t max_results, size_t* n_results) {
    if (!idx || !pattern || pat_len == 0 || !results || !n_results)
        return BAMSI_STATUS_INVALID_ARGUMENT;
    try {
        std::vector<uint8_t> pat(pattern, pattern + pat_len);
        StrandMode mode = static_cast<StrandMode>(idx->idx.header.strand_mode);

        // Forward search
        auto interval = idx->idx.fm.BackwardSearch(pat);
        size_t count = 0;
        for (uint64_t row = interval.lo; row < interval.hi && count < max_results; ++row) {
            if (row == idx->idx.fm.SentinelRow()) continue;
            uint64_t pos = idx->idx.fm.Locate(row);
            if (IsSeparatorPosition(pos, idx->idx.bv.B_read)) continue;

            auto mr = MapOccurrence(pos, pat_len, QueryStrand::Forward,
                                     idx->idx.bv.B_read, idx->idx.reads);
            results[count].chrom_id = mr.chrom_id;
            results[count].p_min = mr.p_min;
            results[count].p_max = mr.p_max;
            results[count].read_id = mr.read_id;
            results[count].is_reverse = 0;
            ++count;
        }

        // Reverse complement search (strand-complete)
        if (mode == StrandMode::StrandComplete) {
            auto rc = ReverseComplement(pat);
            if (rc != pat) {  // avoid double-counting palindromes
                auto rc_interval = idx->idx.fm.BackwardSearch(rc);
                for (uint64_t row = rc_interval.lo;
                     row < rc_interval.hi && count < max_results; ++row) {
                    if (row == idx->idx.fm.SentinelRow()) continue;
                    uint64_t pos = idx->idx.fm.Locate(row);
                    if (IsSeparatorPosition(pos, idx->idx.bv.B_read)) continue;

                    auto mr = MapOccurrence(pos, pat_len, QueryStrand::Reverse,
                                             idx->idx.bv.B_read, idx->idx.reads);
                    results[count].chrom_id = mr.chrom_id;
                    results[count].p_min = mr.p_min;
                    results[count].p_max = mr.p_max;
                    results[count].read_id = mr.read_id;
                    results[count].is_reverse = 1;
                    ++count;
                }
            }
        }

        *n_results = count;
        return BAMSI_STATUS_OK;
    } catch (...) {
        return BAMSI_STATUS_INTERNAL_ERROR;
    }
}

// ─── Query: RegionalCount ──────────────────────────────────────────────

extern "C" bamsi_status_t bamsi_regional_count(const bamsi_index_t* idx,
                                                const uint8_t* pattern, size_t pat_len,
                                                const char* chrom,
                                                uint64_t start, uint64_t end,
                                                uint64_t* count) {
    if (!idx || !pattern || pat_len == 0 || !chrom || !count)
        return BAMSI_STATUS_INVALID_ARGUMENT;
    try {
        std::vector<uint8_t> pat(pattern, pattern + pat_len);
        StrandMode mode = static_cast<StrandMode>(idx->idx.header.strand_mode);
        *count = RegionalCount(pat, std::string(chrom), start, end,
                                idx->idx.fm, idx->idx.bv.B_read, idx->idx.bv.B_window,
                                idx->idx.windows, idx->idx.reads,
                                idx->idx.chrom_to_id, mode);
        return BAMSI_STATUS_OK;
    } catch (...) {
        return BAMSI_STATUS_INTERNAL_ERROR;
    }
}

// ─── Query: RegionalExists ─────────────────────────────────────────────

extern "C" bamsi_status_t bamsi_regional_exists(const bamsi_index_t* idx,
                                                 const uint8_t* pattern, size_t pat_len,
                                                 const char* chrom,
                                                 uint64_t start, uint64_t end,
                                                 uint64_t threshold, int* exists) {
    if (!idx || !pattern || pat_len == 0 || !chrom || !exists)
        return BAMSI_STATUS_INVALID_ARGUMENT;
    try {
        std::vector<uint8_t> pat(pattern, pattern + pat_len);
        StrandMode mode = static_cast<StrandMode>(idx->idx.header.strand_mode);
        *exists = RegionalExists(pat, threshold, std::string(chrom), start, end,
                                  idx->idx.fm, idx->idx.bv.B_read, idx->idx.bv.B_window,
                                  idx->idx.windows, idx->idx.reads,
                                  idx->idx.chrom_to_id, mode) ? 1 : 0;
        return BAMSI_STATUS_OK;
    } catch (...) {
        return BAMSI_STATUS_INTERNAL_ERROR;
    }
}

// ─── Approximate Query Stubs (V1) ──────────────────────────────────────

extern "C" bamsi_status_t bamsi_approx_locate_hamming(const bamsi_index_t* /*idx*/,
                                                       const uint8_t* /*pattern*/,
                                                       size_t /*pat_len*/,
                                                       uint32_t /*max_distance*/,
                                                       bamsi_locate_result_t* /*results*/,
                                                       size_t /*max_results*/,
                                                       size_t* /*n_results*/) {
    return BAMSI_STATUS_NOT_IMPLEMENTED_V1;
}

extern "C" bamsi_status_t bamsi_approx_locate_edit(const bamsi_index_t* /*idx*/,
                                                    const uint8_t* /*pattern*/,
                                                    size_t /*pat_len*/,
                                                    uint32_t /*max_distance*/,
                                                    bamsi_locate_result_t* /*results*/,
                                                    size_t /*max_results*/,
                                                    size_t* /*n_results*/) {
    return BAMSI_STATUS_NOT_IMPLEMENTED_V1;
}
