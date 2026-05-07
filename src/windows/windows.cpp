#include "windows.hpp"

#include <algorithm>
#include <iostream>
#include <map>

namespace bamsi {

uint64_t CigarRefSpan(const CigarRecord& cigar) {
    uint64_t span = 0;
    for (const auto& op : cigar) {
        // M=0, I=1, D=2, N=3, S=4, H=5, P=6, ==7, X=8
        switch (op.op) {
            case 0:  // M
            case 2:  // D
            case 3:  // N
            case 7:  // =
            case 8:  // X
                span += op.len;
                break;
            default:
                break;
        }
    }
    return span;
}

WindowTable BuildWindows(const std::vector<OrderedRead>& reads,
                         const SequenceBundle& bundle,
                         uint64_t T) {
    WindowTable windows;
    const uint64_t N = reads.size();
    if (N == 0) return windows;

    // Group reads by chrom_id (reads are already sorted by chrom_id then pos)
    // Build chromosome → read index ranges
    struct ChromRange {
        uint64_t first_idx;
        uint64_t last_idx;
    };
    std::map<uint32_t, ChromRange> chrom_ranges;

    uint32_t cur_chrom = reads[0].chrom_id;
    uint64_t range_start = 0;
    for (uint64_t i = 0; i < N; ++i) {
        if (reads[i].chrom_id != cur_chrom) {
            chrom_ranges[cur_chrom] = {range_start, i - 1};
            cur_chrom = reads[i].chrom_id;
            range_start = i;
        }
    }
    chrom_ranges[cur_chrom] = {range_start, N - 1};

    // For each chromosome, build windows per Architecture §4.8.2
    for (const auto& [chrom_id, range] : chrom_ranges) {
        uint64_t idx = range.first_idx;
        while (idx <= range.last_idx) {
            uint64_t start_read_idx = idx;
            const auto& start_read = reads[start_read_idx];

            uint64_t l = bundle.readStarts[start_read.read_id];
            uint64_t tentative_end_S = l + T - 1;

            // Extend to include all reads whose S-start falls within [l, tentative_end_S]
            uint64_t last_idx = idx;
            while (last_idx + 1 <= range.last_idx) {
                uint64_t next_start = bundle.readStarts[reads[last_idx + 1].read_id];
                if (next_start <= tentative_end_S) {
                    last_idx++;
                } else {
                    break;
                }
            }

            const auto& last_read = reads[last_idx];

            // r = S-position of the last base of last_read
            uint64_t r = bundle.readStarts[last_read.read_id] +
                         last_read.seq.size() - 1;

            // Include trailing '#' separator if this is not the last read globally
            if (last_read.read_id < N - 1) {
                r = r + 1;  // include the '#'
            }

            // Genomic span
            uint64_t genomic_start = start_read.pos;  // 1-based
            uint64_t genomic_end = 0;
            for (uint64_t j = start_read_idx; j <= last_idx; ++j) {
                uint64_t ref_span = CigarRefSpan(reads[j].cigar);
                uint64_t read_end;
                if (ref_span == 0) {
                    read_end = reads[j].pos;
                } else {
                    read_end = reads[j].pos + ref_span - 1;
                }
                genomic_end = std::max(genomic_end, read_end);
            }

            windows.push_back(Window{
                .chrom_id      = chrom_id,
                .l             = l,
                .r             = r,
                .first_read_id = start_read.read_id,
                .last_read_id  = last_read.read_id,
                .genomic_start = genomic_start,
                .genomic_end   = genomic_end,
            });

            idx = last_idx + 1;
        }
    }

    return windows;
}

}  // namespace bamsi
