/// BAMSIX C ABI Implementation — Contract §10.3
/// All functions are extern "C", thread-safe for concurrent reads.

#include "bamsix/bamsix.h"
#include "bamsix/config.hpp"
#include "format/format.hpp"
#include "query/query.hpp"
#include "mapping/mapping.hpp"

#include <cstring>
#include <memory>
#include <string>
#include <vector>

using namespace bamsix;

/// Opaque index handle — wraps a LoadedIndex.
struct bamsix_index {
    LoadedIndex idx;
};

// ─── Version ───────────────────────────────────────────────────────────

extern "C" const char* bamsix_version(void) {
    return BAMSIX_VERSION;
}

extern "C" uint16_t bamsix_format_version(void) {
    return BAMSIX_FORMAT_VERSION;
}

// ─── Index lifecycle ───────────────────────────────────────────────────

extern "C" bamsix_status_t bamsix_open(const char* path, bamsix_index_t** out) {
    if (!path || !out) return BAMSIX_STATUS_INVALID_ARGUMENT;
    try {
        auto* handle = new bamsix_index();
        handle->idx = ReadBsi(std::string(path));
        *out = handle;
        return BAMSIX_STATUS_OK;
    } catch (const Error& e) {
        // H6 fix: map all ErrorCode variants to C ABI status codes
        switch (e.code) {
            case ErrorCode::CORRUPT_BSI:             return BAMSIX_STATUS_CORRUPT_BSI;
            case ErrorCode::CHECKSUM_MISMATCH:       return BAMSIX_STATUS_CHECKSUM_MISMATCH;
            case ErrorCode::ORDERING_HASH_MISMATCH:  return BAMSIX_STATUS_ORDERING_HASH_MISMATCH;
            case ErrorCode::VERSION_MISMATCH:        return BAMSIX_STATUS_VERSION_MISMATCH;
            case ErrorCode::LOSSY_RECONSTRUCTION:    return BAMSIX_STATUS_LOSSY_RECONSTRUCTION;
            case ErrorCode::UNSUPPORTED_CODEC:       return BAMSIX_STATUS_UNSUPPORTED_CODEC;
            case ErrorCode::STREAM_DECODE_ERROR:     return BAMSIX_STATUS_STREAM_DECODE_ERROR;
            case ErrorCode::INVALID_BAM_INPUT:       return BAMSIX_STATUS_INVALID_ARGUMENT;
            case ErrorCode::BUILD_VALIDATION_FAILED: return BAMSIX_STATUS_BUILD_FAILED;
            default:                                 return BAMSIX_STATUS_INTERNAL_ERROR;
        }
    } catch (...) {
        return BAMSIX_STATUS_INTERNAL_ERROR;
    }
}

extern "C" void bamsix_free(bamsix_index_t** idx) {
    if (idx && *idx) {
        delete *idx;
        *idx = nullptr;
    }
}

// ─── Verification ──────────────────────────────────────────────────────

extern "C" bamsix_status_t bamsix_verify(const char* path, int* valid_out) {
    if (!path || !valid_out) return BAMSIX_STATUS_INVALID_ARGUMENT;
    try {
        *valid_out = VerifyBsi(std::string(path)) ? 1 : 0;
        return BAMSIX_STATUS_OK;
    } catch (...) {
        return BAMSIX_STATUS_INTERNAL_ERROR;
    }
}

// ─── Index info ────────────────────────────────────────────────────────

extern "C" bamsix_status_t bamsix_get_n_reads(const bamsix_index_t* idx, uint64_t* out) {
    if (!idx || !out) return BAMSIX_STATUS_INVALID_ARGUMENT;
    *out = idx->idx.header.N_reads;
    return BAMSIX_STATUS_OK;
}

extern "C" bamsix_status_t bamsix_get_s_length(const bamsix_index_t* idx, uint64_t* out) {
    if (!idx || !out) return BAMSIX_STATUS_INVALID_ARGUMENT;
    *out = idx->idx.header.S_length;
    return BAMSIX_STATUS_OK;
}

extern "C" bamsix_status_t bamsix_get_n_windows(const bamsix_index_t* idx, uint32_t* out) {
    if (!idx || !out) return BAMSIX_STATUS_INVALID_ARGUMENT;
    *out = idx->idx.header.N_windows;
    return BAMSIX_STATUS_OK;
}

extern "C" bamsix_status_t bamsix_get_chrom_count(const bamsix_index_t* idx, uint32_t* out) {
    if (!idx || !out) return BAMSIX_STATUS_INVALID_ARGUMENT;
    *out = idx->idx.header.chrom_count;
    return BAMSIX_STATUS_OK;
}

extern "C" bamsix_status_t bamsix_get_chrom_name(const bamsix_index_t* idx,
                                                uint32_t chrom_idx,
                                                char* buf, size_t buf_len,
                                                size_t* out_len) {
    if (!idx || !out_len) return BAMSIX_STATUS_INVALID_ARGUMENT;
    if (chrom_idx >= idx->idx.chrom_names.size()) return BAMSIX_STATUS_INVALID_ARGUMENT;
    const auto& name = idx->idx.chrom_names[chrom_idx];
    *out_len = name.size();
    if (buf && buf_len > 0) {
        size_t copy_len = std::min(buf_len - 1, name.size());
        std::memcpy(buf, name.c_str(), copy_len);
        buf[copy_len] = '\0';
    }
    return BAMSIX_STATUS_OK;
}

extern "C" bamsix_status_t bamsix_is_lossless(const bamsix_index_t* idx, int* out) {
    if (!idx || !out) return BAMSIX_STATUS_INVALID_ARGUMENT;
    *out = idx->idx.header.is_lossless ? 1 : 0;
    return BAMSIX_STATUS_OK;
}

// ─── Query: GlobalCount ────────────────────────────────────────────────

extern "C" bamsix_status_t bamsix_global_count(const bamsix_index_t* idx,
                                              const uint8_t* pattern, size_t pat_len,
                                              uint64_t* count) {
    if (!idx || !pattern || pat_len == 0 || !count) return BAMSIX_STATUS_INVALID_ARGUMENT;
    try {
        std::vector<uint8_t> pat(pattern, pattern + pat_len);
        StrandMode mode = static_cast<StrandMode>(idx->idx.header.strand_mode);
        *count = GlobalCount(pat, idx->idx.fm, mode);
        return BAMSIX_STATUS_OK;
    } catch (...) {
        return BAMSIX_STATUS_INTERNAL_ERROR;
    }
}

// ─── Query: GlobalExists ───────────────────────────────────────────────

extern "C" bamsix_status_t bamsix_global_exists(const bamsix_index_t* idx,
                                               const uint8_t* pattern, size_t pat_len,
                                               uint64_t threshold, int* exists) {
    if (!idx || !pattern || pat_len == 0 || !exists) return BAMSIX_STATUS_INVALID_ARGUMENT;
    try {
        std::vector<uint8_t> pat(pattern, pattern + pat_len);
        StrandMode mode = static_cast<StrandMode>(idx->idx.header.strand_mode);
        *exists = GlobalExists(pat, idx->idx.fm, mode) ? 1 : 0;
        return BAMSIX_STATUS_OK;
    } catch (...) {
        return BAMSIX_STATUS_INTERNAL_ERROR;
    }
}

// ─── Query: Locate ─────────────────────────────────────────────────────

extern "C" bamsix_status_t bamsix_locate(const bamsix_index_t* idx,
                                        const uint8_t* pattern, size_t pat_len,
                                        bamsix_locate_result_t* results,
                                        size_t max_results, size_t* n_results) {
    if (!idx || !pattern || pat_len == 0 || !results || !n_results)
        return BAMSIX_STATUS_INVALID_ARGUMENT;
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

            auto mr = (!idx->idx.meta_directory.empty() && !idx->idx.map_directory.empty())
                ? MapOccurrenceLazy(pos, pat_len, QueryStrand::Forward,
                                    idx->idx.bv.B_read,
                                    idx->idx.meta_payload, idx->idx.map_payload,
                                    idx->idx.meta_directory, idx->idx.map_directory,
                                    idx->idx.header.meta_codec_id, idx->idx.header.map_codec_id)
                : MapOccurrence(pos, pat_len, QueryStrand::Forward,
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

                    auto mr = (!idx->idx.meta_directory.empty() && !idx->idx.map_directory.empty())
                        ? MapOccurrenceLazy(pos, pat_len, QueryStrand::Reverse,
                                            idx->idx.bv.B_read,
                                            idx->idx.meta_payload, idx->idx.map_payload,
                                            idx->idx.meta_directory, idx->idx.map_directory,
                                            idx->idx.header.meta_codec_id, idx->idx.header.map_codec_id)
                        : MapOccurrence(pos, pat_len, QueryStrand::Reverse,
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
        return BAMSIX_STATUS_OK;
    } catch (...) {
        return BAMSIX_STATUS_INTERNAL_ERROR;
    }
}

// ─── Query: RegionalCount ──────────────────────────────────────────────

extern "C" bamsix_status_t bamsix_regional_count(const bamsix_index_t* idx,
                                                const uint8_t* pattern, size_t pat_len,
                                                const char* chrom,
                                                uint64_t start, uint64_t end,
                                                uint64_t* count) {
    if (!idx || !pattern || pat_len == 0 || !chrom || !count)
        return BAMSIX_STATUS_INVALID_ARGUMENT;
    try {
        std::vector<uint8_t> pat(pattern, pattern + pat_len);
        StrandMode mode = static_cast<StrandMode>(idx->idx.header.strand_mode);
        if (!idx->idx.meta_directory.empty() && !idx->idx.map_directory.empty()) {
            *count = RegionalCountLazy(pat, std::string(chrom), start, end,
                                        idx->idx.fm, idx->idx.bv.B_read, idx->idx.bv.B_window,
                                        idx->idx.windows,
                                        idx->idx.meta_payload, idx->idx.map_payload,
                                        idx->idx.meta_directory, idx->idx.map_directory,
                                        idx->idx.header.meta_codec_id, idx->idx.header.map_codec_id,
                                        idx->idx.chrom_to_id, mode);
        } else {
            *count = RegionalCount(pat, std::string(chrom), start, end,
                                    idx->idx.fm, idx->idx.bv.B_read, idx->idx.bv.B_window,
                                    idx->idx.windows, idx->idx.reads,
                                    idx->idx.chrom_to_id, mode);
        }
        return BAMSIX_STATUS_OK;
    } catch (...) {
        return BAMSIX_STATUS_INTERNAL_ERROR;
    }
}

// ─── Query: RegionalExists ─────────────────────────────────────────────

extern "C" bamsix_status_t bamsix_regional_exists(const bamsix_index_t* idx,
                                                 const uint8_t* pattern, size_t pat_len,
                                                 const char* chrom,
                                                 uint64_t start, uint64_t end,
                                                 uint64_t threshold, int* exists) {
    if (!idx || !pattern || pat_len == 0 || !chrom || !exists)
        return BAMSIX_STATUS_INVALID_ARGUMENT;
    try {
        std::vector<uint8_t> pat(pattern, pattern + pat_len);
        StrandMode mode = static_cast<StrandMode>(idx->idx.header.strand_mode);
        if (!idx->idx.meta_directory.empty() && !idx->idx.map_directory.empty()) {
            *exists = RegionalExistsLazy(pat, threshold, std::string(chrom), start, end,
                                          idx->idx.fm, idx->idx.bv.B_read, idx->idx.bv.B_window,
                                          idx->idx.windows,
                                          idx->idx.meta_payload, idx->idx.map_payload,
                                          idx->idx.meta_directory, idx->idx.map_directory,
                                          idx->idx.header.meta_codec_id, idx->idx.header.map_codec_id,
                                          idx->idx.chrom_to_id, mode) ? 1 : 0;
        } else {
            *exists = RegionalExists(pat, threshold, std::string(chrom), start, end,
                                      idx->idx.fm, idx->idx.bv.B_read, idx->idx.bv.B_window,
                                      idx->idx.windows, idx->idx.reads,
                                      idx->idx.chrom_to_id, mode) ? 1 : 0;
        }
        return BAMSIX_STATUS_OK;
    } catch (...) {
        return BAMSIX_STATUS_INTERNAL_ERROR;
    }
}

// ─── Approximate Query Stubs (V1) ──────────────────────────────────────

extern "C" bamsix_status_t bamsix_approx_locate_hamming(const bamsix_index_t* /*idx*/,
                                                       const uint8_t* /*pattern*/,
                                                       size_t /*pat_len*/,
                                                       uint32_t /*max_distance*/,
                                                       bamsix_locate_result_t* /*results*/,
                                                       size_t /*max_results*/,
                                                       size_t* /*n_results*/) {
    return BAMSIX_STATUS_NOT_IMPLEMENTED_V1;
}

extern "C" bamsix_status_t bamsix_approx_locate_edit(const bamsix_index_t* /*idx*/,
                                                    const uint8_t* /*pattern*/,
                                                    size_t /*pat_len*/,
                                                    uint32_t /*max_distance*/,
                                                    bamsix_locate_result_t* /*results*/,
                                                    size_t /*max_results*/,
                                                    size_t* /*n_results*/) {
    return BAMSIX_STATUS_NOT_IMPLEMENTED_V1;
}

// ─── Streaming Locate Iterator (Contract §10.3) ────────────────────────

struct bamsix_locate_iter {
    std::vector<Match> matches;
    std::vector<std::string> chrom_names;
    size_t cursor = 0;
};

extern "C" bamsix_status_t bamsix_locate_iter_create(const bamsix_index_t* idx,
                                                    const uint8_t* pattern,
                                                    size_t pat_len,
                                                    bamsix_locate_iter_t** iter) {
    if (!idx || !pattern || !iter) return BAMSIX_STATUS_INVALID_ARGUMENT;
    if (pat_len == 0) return BAMSIX_STATUS_INVALID_ARGUMENT;
    try {
        std::vector<uint8_t> P(pattern, pattern + pat_len);
        StrandMode mode = static_cast<StrandMode>(idx->idx.header.strand_mode);

        auto matches = Locate(P, idx->idx.fm, idx->idx.bv.B_read,
                              idx->idx.reads, idx->idx.chrom_names,
                              mode, false /* unsorted for streaming */);

        auto* it = new bamsix_locate_iter();
        it->matches = std::move(matches);
        it->chrom_names = idx->idx.chrom_names;
        it->cursor = 0;
        *iter = it;
        return BAMSIX_STATUS_OK;
    } catch (const Error& e) {
        if (e.code == ErrorCode::INVALID_PATTERN) return BAMSIX_STATUS_INVALID_ARGUMENT;
        if (e.code == ErrorCode::EMPTY_PATTERN) return BAMSIX_STATUS_INVALID_ARGUMENT;
        return BAMSIX_STATUS_INTERNAL_ERROR;
    } catch (...) {
        return BAMSIX_STATUS_INTERNAL_ERROR;
    }
}

extern "C" bamsix_status_t bamsix_locate_iter_next(bamsix_locate_iter_t* iter,
                                                  bamsix_locate_result_t* result,
                                                  int* has_more) {
    if (!iter || !result || !has_more) return BAMSIX_STATUS_INVALID_ARGUMENT;

    if (iter->cursor >= iter->matches.size()) {
        *has_more = 0;
        return BAMSIX_STATUS_OK;
    }

    const auto& m = iter->matches[iter->cursor];
    // Find chrom_id from name
    result->chrom_id = 0;
    for (uint32_t ci = 0; ci < iter->chrom_names.size(); ++ci) {
        if (iter->chrom_names[ci] == m.chrom) {
            result->chrom_id = ci;
            break;
        }
    }
    result->p_min = m.p_min;
    result->p_max = m.p_max;
    result->read_id = m.read_id;
    result->is_reverse = (m.query_strand == QueryStrand::Reverse) ? 1 : 0;

    iter->cursor++;
    *has_more = (iter->cursor < iter->matches.size()) ? 1 : 0;
    return BAMSIX_STATUS_OK;
}

extern "C" void bamsix_locate_iter_free(bamsix_locate_iter_t** iter) {
    if (iter && *iter) {
        delete *iter;
        *iter = nullptr;
    }
}
