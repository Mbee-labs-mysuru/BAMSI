#include "../src/format/format.hpp"
#include "../src/query/query.hpp"
#include <iostream>
using namespace bamsix;
int main() {
    auto idx = ReadBsi("data/test/synthetic_10reads.bsi");
    std::cerr << "Loaded: |S|=" << idx.header.S_length
              << " reads=" << idx.header.N_reads
              << " sentinel_row=" << idx.header.sentinel_row << "\n";
    
    // Check FM basics
    std::cerr << "FM: bwt_len=" << idx.fm.SLen()
              << " sentinel=" << idx.fm.SentinelRow()
              << " SA_samples=" << idx.fm.SASamples().size() << "\n";
    
    // Print SA samples
    std::cerr << "SA_samples:";
    for (auto s : idx.fm.SASamples()) std::cerr << " " << s;
    std::cerr << "\n";
    
    // Test count
    std::vector<uint8_t> pat = {0,1,2,3}; // ACGT
    auto interval = idx.fm.BackwardSearch(pat);
    std::cerr << "ACGT interval: [" << interval.lo << ", " << interval.hi << ") = " << interval.size() << "\n";
    
    // Test locate for first few
    uint64_t n1 = idx.header.S_length + 1; // |S$| = |S|+1
    for (uint64_t row = interval.lo; row < interval.hi && row < interval.lo + 5; row++) {
        if (row == idx.fm.SentinelRow()) continue;
        uint64_t pos_raw = idx.fm.Locate(row);
        uint64_t pos = pos_raw % n1;
        std::cerr << "  row " << row << " raw_pos=" << pos_raw << " mod_pos=" << pos << "\n";
    }
    return 0;
}
