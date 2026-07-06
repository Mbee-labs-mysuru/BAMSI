#pragma once

#include "bamsix/types.hpp"
#include "../fmindex/fmindex.hpp"
#include "../bitvectors/bitvectors.hpp"
#include "../streamencode/streamencode.hpp"
#include "../sarange/sarange.hpp"
#include <string>

namespace bamsix {

/// All data required to write a sealed .bsi file.
struct SealInput {
    BsiHeader                    header;
    SeqEncodeResult              seq;
    QualEncodeResult             qual;
    MetaEncodeResult             meta;
    MapEncodeResult              map;
    FMIndexEngine                fm;
    BitVectors                   bv;
    WindowTable                  windows;
    std::vector<OrderedRead>     reads;  // for query-time Locate

    // V5 ENHANCED features
    std::vector<uint64_t>        isa_samples;   // ISA samples (§4.4)
    uint64_t                     isa_step = 0;  // ISA sample step s'
    SARange                      sarange;       // SARange wavelet tree (§5.3)
};

/// Write a sealed .bsi file per Architecture §7.
/// Writes header, stream sections, FM-index section, bitvector section,
/// window section, directory section, ISA/SARange sections, and footer.
void WriteBsi(const std::string& path, const SealInput& input);

}  // namespace bamsix

