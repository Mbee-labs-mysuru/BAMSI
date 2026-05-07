#pragma once

#include "bamsi/types.hpp"
#include <vector>

namespace bamsi {

/// Result of SA-IS construction (Stage 4 of build pipeline).
struct SaisResult {
    std::vector<int64_t>  SA;           // suffix array
    std::vector<uint8_t>  BWT;          // BWT[i] = S[(SA[i]-1+|S|) mod |S|]
    uint64_t              sentinel_row; // row index whose SA value would be |S| (the $)
    std::vector<uint64_t> ISA_samples;  // ISA_samples[k] = ISA[k * s'] (Architecture §4.4)
};

/// Run SA-IS on the sequence S to produce the suffix array and BWT.
/// The sentinel '$' is conceptual (not stored in S); SA-IS uses the integer
/// alphabet {0..5} directly (codes A,C,G,T,N,#).
/// Per Architecture §4.4: SA-IS uses the Nong, Zhang & Chan (2009) algorithm.
SaisResult ComputeSuffixArray(const SequenceBundle& bundle);

/// Compute ISA samples from SA (Architecture §4.4 step 5).
/// ISA[SA[j]] = j for all j, then sample ISA_samples[k] = ISA[k * s'].
/// Full ISA is not stored; only the sampled array.
void ComputeISASamples(SaisResult& result, uint64_t sample_step_s_prime);

}  // namespace bamsi
