#include "sais.hpp"

#include <libsais.h>
#include <libsais64.h>
#include <cstring>
#include <stdexcept>

namespace bamsi {

SaisResult ComputeSuffixArray(const SequenceBundle& bundle) {
    SaisResult result;
    const auto& S = bundle.S;
    const int64_t n = static_cast<int64_t>(S.size());

    if (n == 0) {
        return result;
    }

    // ─── Stage 4 (Architecture §4.4) ────────────────────────────────────────
    // Conceptually form S$ with |S|+1 characters.
    // $ is lexicographically smaller than all codes 0–5.
    //
    // Implementation: shift all codes up by 1 so sentinel = 0.
    // A=1, C=2, G=3, T=4, N=5, #=6.
    std::vector<uint8_t> input(n + 1);
    for (int64_t i = 0; i < n; ++i) {
        input[i] = S[i] + 1;  // shift: 0→1, 1→2, ..., 5→6
    }
    input[n] = 0;  // sentinel (smallest symbol)

    int64_t sa_len = n + 1;

    if (sa_len <= INT32_MAX) {
        std::vector<int32_t> sa32(sa_len);
        // libsais(T, SA, n, fs, freq): fs=0, no extra space
        int ret = libsais(input.data(), sa32.data(), static_cast<int32_t>(sa_len),
                          0, nullptr);
        if (ret != 0) {
            throw Error{ErrorCode::BUILD_VALIDATION_FAILED,
                        "libsais SA-IS construction failed"};
        }

        result.SA.resize(sa_len);
        result.sentinel_row = 0;
        for (int64_t i = 0; i < sa_len; ++i) {
            result.SA[i] = sa32[i];
            if (sa32[i] == static_cast<int32_t>(n)) {
                result.sentinel_row = static_cast<uint64_t>(i);
            }
        }
    } else {
        std::vector<int64_t> sa64(sa_len);
        int ret = libsais64(input.data(), sa64.data(), sa_len, 0, nullptr);
        if (ret != 0) {
            throw Error{ErrorCode::BUILD_VALIDATION_FAILED,
                        "libsais64 SA-IS construction failed"};
        }
        result.SA = std::move(sa64);
        for (int64_t i = 0; i < sa_len; ++i) {
            if (result.SA[i] == n) {
                result.sentinel_row = static_cast<uint64_t>(i);
                break;
            }
        }
    }

    // ─── Derive BWT (Architecture §4.4) ─────────────────────────────────────
    // BWT[i] = S[(SA[i] - 1 + n) mod n]   for i != sentinel_row
    // BWT[sentinel_row] = CODE_SENT (conceptual $, code 6)
    //
    // We store ALL |S|+1 entries in the BWT, with the sentinel row
    // explicitly marked as CODE_SENT. The FM-index handles this internally.
    // Architecture §4.4 says "$ is not stored" in the .bsi file, but
    // in-memory during build, we keep it for correct FM construction.
    // The seal module will strip it when writing to disk.
    result.BWT.resize(static_cast<size_t>(sa_len));

    for (int64_t i = 0; i < sa_len; ++i) {
        if (static_cast<uint64_t>(i) == result.sentinel_row) {
            // BWT[sentinel_row] = S$[(SA[sentinel]-1+|S$|) mod |S$|]
            //                   = S$[(n-1+n+1) mod (n+1)]
            //                   = S$[n-1] = S[n-1]  (last char of S)
            result.BWT[i] = S[n - 1];
            continue;
        }

        int64_t sa_val = result.SA[i];
        if (sa_val == 0) {
            // BWT[i] = S$[(0-1+n+1) mod (n+1)] = S$[n] = $ (sentinel)
            result.BWT[i] = CODE_SENT;
        } else {
            result.BWT[i] = S[sa_val - 1];
        }
    }

    return result;
}

void ComputeISASamples(SaisResult& result, uint64_t sample_step_s_prime) {
    if (sample_step_s_prime == 0) return;

    const int64_t sa_len = static_cast<int64_t>(result.SA.size());
    if (sa_len == 0) return;

    // Step 1: Compute full ISA by inverting SA: ISA[SA[j]] = j
    // We only need ISA for text positions [0..|S|-1], not the sentinel position |S|.
    const int64_t n = sa_len - 1;  // |S| (SA has |S|+1 entries for S$)
    std::vector<uint64_t> isa(n);

    for (int64_t j = 0; j < sa_len; ++j) {
        int64_t sa_val = result.SA[j];
        if (sa_val < n) {  // skip sentinel position (SA[j] == |S|)
            isa[sa_val] = static_cast<uint64_t>(j);
        }
    }

    // Step 2: Sample ISA at intervals of s'
    // ISA_samples[k] = ISA[k * s'] for k = 0..floor(|S| / s')
    uint64_t n_samples = static_cast<uint64_t>(n) / sample_step_s_prime + 1;
    result.ISA_samples.resize(n_samples);

    for (uint64_t k = 0; k < n_samples; ++k) {
        uint64_t text_pos = k * sample_step_s_prime;
        if (text_pos < static_cast<uint64_t>(n)) {
            result.ISA_samples[k] = isa[text_pos];
        } else {
            result.ISA_samples[k] = 0;
        }
    }
    // Full ISA vector is now discarded (goes out of scope)
}

}  // namespace bamsi
