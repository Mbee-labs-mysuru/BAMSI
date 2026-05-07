#include "seqbuilder.hpp"

#include <numeric>

namespace bamsi {

SequenceBundle BuildSequence(const std::vector<OrderedRead>& reads) {
    SequenceBundle bundle;
    const uint64_t N = reads.size();
    if (N == 0) {
        return bundle;
    }

    // Compute total length: |S| = sum |r_i| + (N - 1) separators
    uint64_t total_len = 0;
    for (const auto& r : reads) {
        total_len += r.seq.size();
    }
    total_len += (N - 1);  // N-1 '#' separators

    bundle.S.reserve(total_len);
    bundle.readStarts.resize(N);

    for (uint64_t i = 0; i < N; ++i) {
        bundle.readStarts[i] = bundle.S.size();

        // Append read bases
        bundle.S.insert(bundle.S.end(), reads[i].seq.begin(), reads[i].seq.end());

        // Append '#' separator after every read except the last
        if (i + 1 < N) {
            bundle.S.push_back(CODE_SEP);
        }
    }

    // Verify length
    if (bundle.S.size() != total_len) {
        throw Error{ErrorCode::BUILD_VALIDATION_FAILED,
                    "Sequence length mismatch: expected " + std::to_string(total_len) +
                    " got " + std::to_string(bundle.S.size())};
    }

    return bundle;
}

}  // namespace bamsi
