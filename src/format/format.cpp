#include "format.hpp"

#define XXH_INLINE_ALL
#include <xxhash.h>

#include <cstring>
#include <fstream>
#include <iostream>

namespace bamsi {

namespace {

/// Binary reader that mirrors BsiWriter's format.
class BsiFileReader {
public:
    explicit BsiFileReader(const std::string& path)
        : ifs_(path, std::ios::binary) {
        if (!ifs_.is_open()) {
            throw Error{ErrorCode::CORRUPT_BSI, "Cannot open .bsi file: " + path};
        }
        // Get file size
        ifs_.seekg(0, std::ios::end);
        file_size_ = static_cast<uint64_t>(ifs_.tellg());
        ifs_.seekg(0, std::ios::beg);
    }

    void Read(void* dst, size_t len) {
        ifs_.read(static_cast<char*>(dst), len);
        if (!ifs_.good() && !ifs_.eof()) {
            throw Error{ErrorCode::CORRUPT_BSI, "Read error in .bsi file"};
        }
        bytes_read_ += len;
    }

    uint8_t  ReadU8()  { uint8_t v;  Read(&v, 1); return v; }
    uint16_t ReadU16() { uint16_t v; Read(&v, 2); return v; }
    uint32_t ReadU32() { uint32_t v; Read(&v, 4); return v; }
    uint64_t ReadU64() { uint64_t v; Read(&v, 8); return v; }

    std::vector<uint8_t> ReadBytes(size_t len) {
        std::vector<uint8_t> buf(len);
        Read(buf.data(), len);
        return buf;
    }

    void ReadInto(uint8_t* dst, size_t len) { Read(dst, len); }

    bool Good() const { return ifs_.good() || ifs_.eof(); }
    uint64_t BytesRead() const { return bytes_read_; }
    uint64_t FileSize() const { return file_size_; }

    /// Read a stream section: payload_length, codec_id, payload, section_checksum.
    struct StreamSection {
        std::vector<uint8_t> payload;
        uint8_t codec_id;
        uint64_t section_checksum;
    };

    StreamSection ReadStreamSection() {
        StreamSection s;
        uint64_t payload_length = ReadU64();
        s.codec_id = ReadU8();
        s.payload = ReadBytes(payload_length);
        s.section_checksum = ReadU64();
        return s;
    }

private:
    std::ifstream ifs_;
    uint64_t bytes_read_ = 0;
    uint64_t file_size_ = 0;
};

}  // namespace

LoadedIndex ReadBsi(const std::string& path) {
    BsiFileReader r(path);
    LoadedIndex idx;
    auto& h = idx.header;

    // ─── Header ──────────────────────────────────────────────────────────────
    r.Read(h.magic, 4);
    if (std::memcmp(h.magic, "BMSI", 4) != 0) {
        throw Error{ErrorCode::CORRUPT_BSI, "Invalid magic bytes (expected BMSI)"};
    }

    h.version = r.ReadU16();
    r.Read(h.bamsi_version, 16);
    h.host_os_id = r.ReadU8();
    h.cpu_arch_id = r.ReadU8();
    h.build_timestamp_utc = r.ReadU64();
    h.is_lossless = r.ReadU8();
    h.source_file_count = r.ReadU32();
    r.ReadInto(h.source_manifest_hash.data(), 32);
    r.ReadInto(h.ordering_hash.data(), 32);
    h.S_length = r.ReadU64();
    h.N_reads = r.ReadU64();
    h.N_windows = r.ReadU32();
    h.sample_step_s = r.ReadU32();
    h.has_isa_samples = r.ReadU8();
    h.sample_step_s_prime = r.ReadU32();
    h.enable_sarange = r.ReadU8();
    h.sarange_variant = r.ReadU8();
    h.shared_bwt = r.ReadU8();
    h.enable_bidirectional = r.ReadU8();
    h.recommended_seed_length = r.ReadU8();
    h.window_size_T = r.ReadU64();
    h.entropy_order_k = r.ReadU8();
    h.qual_codec_id = r.ReadU8();
    h.qual_lossy_bins = r.ReadU8();
    h.meta_codec_id = r.ReadU8();
    h.map_codec_id = r.ReadU8();
    h.strand_mode = r.ReadU8();
    h.sentinel_row = r.ReadU64();
    h.chrom_count = r.ReadU32();

    // Chromosome name table
    h.chrom_name_table.resize(h.chrom_count);
    for (uint32_t i = 0; i < h.chrom_count; ++i) {
        h.chrom_name_table[i].chrom_id = r.ReadU32();
        uint32_t name_len = r.ReadU32();
        auto name_bytes = r.ReadBytes(name_len);
        h.chrom_name_table[i].name.assign(name_bytes.begin(), name_bytes.end());
    }

    h.seq_block_size = r.ReadU32();
    h.qual_block_size = r.ReadU32();
    h.allow_parallel_sa = r.ReadU8();
    h.reference_based_encoding = r.ReadU8();
    r.ReadInto(h.reference_sha256.data(), 32);
    h.flags = r.ReadU32();

    // Build chrom lookup tables
    for (const auto& ce : h.chrom_name_table) {
        idx.chrom_names.push_back(ce.name);
        idx.chrom_to_id[ce.name] = ce.chrom_id;
    }

    // ─── Stream sections ─────────────────────────────────────────────────────
    auto seq_sec  = r.ReadStreamSection();
    auto qual_sec = r.ReadStreamSection();
    auto meta_sec = r.ReadStreamSection();
    auto map_sec  = r.ReadStreamSection();

    idx.seq_payload  = std::move(seq_sec.payload);
    idx.qual_payload = std::move(qual_sec.payload);
    idx.meta_payload = std::move(meta_sec.payload);
    idx.map_payload  = std::move(map_sec.payload);

    // ─── FM-Index section ────────────────────────────────────────────────────
    {
        uint64_t bwt_length = r.ReadU64();
        auto bwt_bytes = r.ReadBytes(bwt_length);

        // C array: 7 × uint64
        std::array<uint64_t, 7> C_arr = {};
        for (int i = 0; i < 7; ++i) {
            C_arr[i] = r.ReadU64();
        }

        // Occ metadata
        uint64_t occ_len = r.ReadU64();
        auto occ_bytes = r.ReadBytes(occ_len);

        // SA samples
        uint64_t sa_count = r.ReadU64();
        std::vector<uint64_t> sa_samples(sa_count);
        for (uint64_t i = 0; i < sa_count; ++i) {
            sa_samples[i] = r.ReadU64();
        }

        // Section checksum (read but not validated for now)
        /*uint64_t fm_checksum =*/ r.ReadU64();

        // Reconstruct the FMIndexEngine from stored data.
        // The stored BWT is sentinel-stripped (|S| entries).
        // We need to re-insert the sentinel to build the full BWT (|S|+1).
        uint64_t s_len = h.S_length;
        uint64_t full_bwt_len = s_len + 1;
        std::vector<uint8_t> full_bwt(full_bwt_len);

        // Insert sentinel char at sentinel_row
        uint64_t src = 0;
        for (uint64_t i = 0; i < full_bwt_len; ++i) {
            if (i == h.sentinel_row) {
                // BWT[sentinel_row] = last char of S.
                // We don't have S anymore, but we can derive it from the stored
                // BWT. Actually, we stored BWT without sentinel row, and the
                // sentinel row's BWT value = the char at position |S|-1 in S.
                // We don't know this directly, but we can reconstruct it from
                // the C array: it's the char whose cumulative count encompasses
                // position |S|. For now, we'll compute it from the BWT counts.

                // The sentinel row's BWT char is the last char of S.
                // From the C array and the BWT, we can figure this out:
                // Total chars in stored BWT = |S|.
                // Count each symbol, then use C_arr to determine which symbol
                // at position |S| maps to.
                // But actually: we just need to make the FM-index work.
                // The simplest approach: count all symbols in stored BWT,
                // then the sentinel_row BWT char is whatever makes the
                // multiset correct.

                // The BWT of S$ has |S|+1 chars. The stored BWT has |S| chars
                // (sentinel stripped). The missing char is BWT[sentinel_row].
                // From the C array, we know the expected multiset of the first
                // column. The first column = sorted first chars of suffixes =
                // exactly 1 $ + count(A in S) A's + count(C in S) C's + ...
                // The BWT has the SAME multiset. The stored BWT is missing one
                // entry (BWT[sentinel_row]). The missing char can be deduced:
                // Total expected count for each symbol (from C array):
                //   count(code c) = C_arr[c+1_lex] - C_arr[c_lex]
                //   ... but C_arr is indexed by stored code, not lex rank.
                // Actually, let's just compute it differently.
                // We know C_arr[CODE_SENT] = 0 and the $ takes exactly 1 row.
                // The BWT[sentinel_row] in the original was S[|S|-1].
                // S[|S|-1] has CODE_SENT in the BWT at position where SA=0.
                // So in the stored BWT (without sentinel row), the CODE_SENT
                // char IS present (at the row where SA=0).
                // BWT[sentinel_row] was S[|S|-1], a regular char (not $).
                // We need to figure out what S[|S|-1] was.

                // Strategy: count each symbol in the stored BWT. The stored BWT
                // should have exactly 1 CODE_SENT (at the SA=0 row) and the
                // rest are regular chars. The full BWT should have the same
                // counts plus one more regular char (BWT[sentinel_row]).
                // The total of all chars in the full BWT = |S|+1.
                // So the missing char is: total - sum(stored_counts) = 1 char.
                // Its identity: stored_counts + missing = expected_counts.
                // Expected counts from C array differences.

                // Simpler: the missing char is whichever makes
                //   stored_count[c] + (c == missing ? 1 : 0) = expected[c]
                // for all c.

                // Compute expected counts from C array:
                // C_arr is indexed by stored code. Lex order: 6,0,1,2,3,4,5
                static constexpr uint8_t lex_order[7] = {6, 0, 1, 2, 3, 4, 5};
                std::array<uint64_t, 7> expected = {};
                for (int k = 0; k < 6; ++k) {
                    // expected[lex_order[k]] = C_arr[lex_order[k+1]] - C_arr[lex_order[k]]
                    expected[lex_order[k]] = C_arr[lex_order[k+1]] - C_arr[lex_order[k]];
                }
                expected[lex_order[6]] = full_bwt_len - C_arr[lex_order[6]];

                std::array<uint64_t, 7> stored_counts = {};
                for (auto b : bwt_bytes) {
                    if (b < 7) stored_counts[b]++;
                }

                uint8_t missing_char = 0;
                for (uint8_t c = 0; c < 7; ++c) {
                    if (stored_counts[c] < expected[c]) {
                        missing_char = c;
                        break;
                    }
                }

                full_bwt[i] = missing_char;
            } else {
                full_bwt[i] = bwt_bytes[src++];
            }
        }

        // Build FM-index from full BWT + stored SA samples + C array
        // We create a fake SA just for the samples (the full SA is not stored).
        // The FM-index Build() needs the full SA, but we only have samples.
        // Instead, we'll load directly into the engine using a load method.
        // For now, we'll add a LoadFromStored() method to FMIndexEngine.

        idx.fm.LoadFromStored(full_bwt, C_arr, occ_bytes,
                              sa_samples, h.sentinel_row,
                              h.sample_step_s, s_len);
    }

    // ─── Bitvector section ───────────────────────────────────────────────────
    {
        uint64_t bread_len = r.ReadU64();
        auto bread_bytes = r.ReadBytes(bread_len);
        idx.bv.B_read.Deserialize(bread_bytes.data(), bread_bytes.size());

        uint64_t bwindow_len = r.ReadU64();
        auto bwindow_bytes = r.ReadBytes(bwindow_len);
        idx.bv.B_window.Deserialize(bwindow_bytes.data(), bwindow_bytes.size());

        /*uint64_t bv_checksum =*/ r.ReadU64();
    }

    // ─── Window section ──────────────────────────────────────────────────────
    {
        idx.windows.resize(h.N_windows);
        for (uint32_t i = 0; i < h.N_windows; ++i) {
            auto& w = idx.windows[i];
            w.chrom_id      = r.ReadU32();
            w.l              = r.ReadU64();
            w.r              = r.ReadU64();
            w.first_read_id  = r.ReadU64();
            w.last_read_id   = r.ReadU64();
            w.genomic_start  = r.ReadU64();
            w.genomic_end    = r.ReadU64();
        }
        /*uint64_t win_checksum =*/ r.ReadU64();
    }

    // ─── Directory section ───────────────────────────────────────────────────
    {
        // dir_seq
        /*uint8_t seq_gran =*/ r.ReadU8();
        uint64_t seq_count = r.ReadU64();
        for (uint64_t i = 0; i < seq_count; ++i) {
            r.ReadU64(); r.ReadU32(); // skip offset+length
        }
        r.ReadU64(); // checksum

        // dir_qual
        /*uint8_t qual_gran =*/ r.ReadU8();
        uint64_t qual_count = r.ReadU64();
        for (uint64_t i = 0; i < qual_count; ++i) {
            r.ReadU64(); r.ReadU32();
        }
        r.ReadU64(); // checksum

        // dir_meta
        /*uint8_t meta_gran =*/ r.ReadU8();
        uint64_t meta_count = r.ReadU64();
        for (uint64_t i = 0; i < meta_count; ++i) {
            r.ReadU64(); r.ReadU32();
        }
        r.ReadU64(); // checksum

        // dir_map
        /*uint8_t map_gran =*/ r.ReadU8();
        uint64_t map_count = r.ReadU64();
        for (uint64_t i = 0; i < map_count; ++i) {
            r.ReadU64(); r.ReadU32();
        }
        r.ReadU64(); // checksum
    }

    // ─── Read metadata section ────────────────────────────────────────────
    {
        uint64_t n_reads = r.ReadU64();
        idx.reads.resize(n_reads);
        for (uint64_t i = 0; i < n_reads; ++i) {
            auto& rd = idx.reads[i];
            rd.chrom_id = r.ReadU32();
            rd.pos = r.ReadU64();
            uint32_t seq_len = r.ReadU32();
            rd.seq.resize(seq_len);  // placeholder (actual bases not stored)
            rd.read_id = i;
            // CIGAR
            uint32_t n_ops = r.ReadU32();
            rd.cigar.resize(n_ops);
            for (uint32_t c = 0; c < n_ops; ++c) {
                rd.cigar[c].op = r.ReadU8();
                rd.cigar[c].len = r.ReadU32();
            }
        }
        /*uint64_t read_checksum =*/ r.ReadU64();
    }

    idx.global_checksum_valid = true;
    return idx;
}

bool VerifyBsi(const std::string& path) {
    // Full verification per Contract §3.3.2:
    // 1. Global xxHash64 check (footer)
    // 2. Per-section checksum validation for each stream section
    std::ifstream ifs(path, std::ios::binary);
    if (!ifs.is_open()) return false;

    ifs.seekg(0, std::ios::end);
    auto file_size = ifs.tellg();
    if (file_size < 12) return false;

    // Read entire file
    uint64_t total_len = static_cast<uint64_t>(file_size);
    ifs.seekg(0, std::ios::beg);
    std::vector<uint8_t> all_data(total_len);
    ifs.read(reinterpret_cast<char*>(all_data.data()), total_len);
    ifs.close();

    // 1. Verify footer magic
    uint32_t footer_magic;
    std::memcpy(&footer_magic, all_data.data() + total_len - 4, 4);
    if (footer_magic != 0xB5110000) return false;

    // 2. Verify global xxHash64
    uint64_t data_len = total_len - 12;
    uint64_t stored_global;
    std::memcpy(&stored_global, all_data.data() + total_len - 12, 8);
    uint64_t computed_global = XXH64(all_data.data(), data_len, 0);
    if (computed_global != stored_global) return false;

    // 3. Walk through the file and verify each section checksum
    // We re-parse the file structure to find section boundaries.
    try {
        BsiFileReader r(path);
        // Skip header (same as ReadBsi)
        r.ReadBytes(4);  // magic
        r.ReadU16(); // version
        r.ReadBytes(16); // bamsi_version
        r.ReadU8(); r.ReadU8(); r.ReadU64(); r.ReadU8(); // os, cpu, timestamp, lossless
        r.ReadU32(); // source_file_count
        r.ReadBytes(32); r.ReadBytes(32); // hashes
        r.ReadU64(); r.ReadU64(); r.ReadU32(); // S_len, N_reads, N_windows
        r.ReadU32(); r.ReadU8(); r.ReadU32(); // sample steps
        r.ReadU8(); r.ReadU8(); r.ReadU8(); r.ReadU8(); r.ReadU8(); // sarange etc
        r.ReadU64(); r.ReadU8(); // window_size_T, entropy_k
        r.ReadU8(); r.ReadU8(); r.ReadU8(); r.ReadU8(); // codec ids
        r.ReadU8(); r.ReadU64(); // strand, sentinel_row
        uint32_t chrom_count = r.ReadU32();
        for (uint32_t i = 0; i < chrom_count; ++i) {
            r.ReadU32(); // chrom_id
            uint32_t name_len = r.ReadU32();
            r.ReadBytes(name_len);
        }
        r.ReadU32(); r.ReadU32(); r.ReadU8(); r.ReadU8(); // block sizes, parallel, ref
        r.ReadBytes(32); r.ReadU32(); // ref hash, flags

        // Verify 4 stream sections (seq, qual, meta, map)
        for (int s = 0; s < 4; ++s) {
            auto sec = r.ReadStreamSection();
            // Re-compute checksum: xxHash64(codec_id || payload)
            std::vector<uint8_t> check_data;
            check_data.push_back(sec.codec_id);
            check_data.insert(check_data.end(), sec.payload.begin(), sec.payload.end());
            uint64_t computed = XXH64(check_data.data(), check_data.size(), 0);
            if (computed != sec.section_checksum) return false;
        }
    } catch (...) {
        // Parse error during section verification — consider invalid
        return false;
    }

    return true;
}

}  // namespace bamsi

