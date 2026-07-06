#include "sais.hpp"

#include <libsais.h>
#include <libsais64.h>
#include <cstring>
#include <stdexcept>

namespace bamsix {

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

    const int64_t n = sa_len - 1;  // |S| (SA has |S|+1 entries for S$)

    // Step 1: Allocate ISA_samples
    // ISA_samples[k] = ISA[k * s'] for k = 0..floor(|S| / s')
    uint64_t n_samples = static_cast<uint64_t>(n) / sample_step_s_prime + 1;
    result.ISA_samples.assign(n_samples, 0);

    // Step 2: Iterate SA, and when SA[j] is a multiple of s', store j.
    // This avoids allocating a full 16GB ISA array for 2B reads.
    for (int64_t j = 0; j < sa_len; ++j) {
        int64_t sa_val = result.SA[j];
        if (sa_val < n && (sa_val % sample_step_s_prime) == 0) {
            result.ISA_samples[sa_val / sample_step_s_prime] = static_cast<uint64_t>(j);
        }
    }
}

}  // namespace bamsix
