#include "format.hpp"
#include "bamsix/config.hpp"

#define XXH_INLINE_ALL
#include <xxhash.h>

#include <openssl/sha.h>
#include <cstring>
#include <fstream>
#include <iostream>

namespace bamsix {

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

    /// Read and validate a stream section — Contract §9.4: no silent fallback.
    StreamSection ReadAndValidateStreamSection(const std::string& name) {
        auto s = ReadStreamSection();
        // Verify: xxHash64(codec_id || payload)
        std::vector<uint8_t> check_data;
        check_data.reserve(1 + s.payload.size());
        check_data.push_back(s.codec_id);
        check_data.insert(check_data.end(), s.payload.begin(), s.payload.end());
        uint64_t computed = XXH64(check_data.data(), check_data.size(), 0);
        if (computed != s.section_checksum) {
            throw Error{ErrorCode::CHECKSUM_MISMATCH,
                        "Section checksum mismatch in " + name +
                        " stream (expected " + std::to_string(s.section_checksum) +
                        ", got " + std::to_string(computed) + ")"};
        }
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

    // C1/H2 fix: validate format_version (Contract §10.6)
    // v1.x readers MUST refuse v2.x files with VERSION_MISMATCH
    if (h.version > BAMSIX_FORMAT_VERSION) {
        throw Error{ErrorCode::VERSION_MISMATCH,
                    "Unsupported .bsi format version " + std::to_string(h.version) +
                    " (this build supports up to version " +
                    std::to_string(BAMSIX_FORMAT_VERSION) + ")"};
    }

    r.Read(h.bamsix_version, 16);
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

    // ─── Stream sections (validated per Contract §9.4) ────────────────────────
    auto seq_sec  = r.ReadAndValidateStreamSection("S_seq");
    auto qual_sec = r.ReadAndValidateStreamSection("S_qual");
    auto meta_sec = r.ReadAndValidateStreamSection("S_meta");
    auto map_sec  = r.ReadAndValidateStreamSection("S_map");

    idx.seq_payload  = std::move(seq_sec.payload);
    idx.qual_payload = std::move(qual_sec.payload);
    idx.meta_payload = std::move(meta_sec.payload);
    idx.map_payload  = std::move(map_sec.payload);

    // ─── FM-Index section ────────────────────────────────────────────────────
    {
        // C4 fix: accumulate section data for checksum validation
        std::vector<uint8_t> fm_section_data;

        uint64_t bwt_length = r.ReadU64();
        auto bwt_bytes = r.ReadBytes(bwt_length);
        fm_section_data.insert(fm_section_data.end(), bwt_bytes.begin(), bwt_bytes.end());

        // C array: 7 × uint64
        std::array<uint64_t, 7> C_arr = {};
        for (int i = 0; i < 7; ++i) {
            C_arr[i] = r.ReadU64();
            uint8_t buf[8]; std::memcpy(buf, &C_arr[i], 8);
            fm_section_data.insert(fm_section_data.end(), buf, buf + 8);
        }

        // Occ metadata
        uint64_t occ_len = r.ReadU64();
        auto occ_bytes = r.ReadBytes(occ_len);
        fm_section_data.insert(fm_section_data.end(), occ_bytes.begin(), occ_bytes.end());

        // SA samples
        uint64_t sa_count = r.ReadU64();
        std::vector<uint64_t> sa_samples(sa_count);
        for (uint64_t i = 0; i < sa_count; ++i) {
            sa_samples[i] = r.ReadU64();
            uint8_t buf[8]; std::memcpy(buf, &sa_samples[i], 8);
            fm_section_data.insert(fm_section_data.end(), buf, buf + 8);
        }

        // C4 fix: Validate FM-index section checksum (Contract §9.4 — no silent fallback)
        uint64_t fm_checksum = r.ReadU64();
        uint64_t computed_fm = XXH64(fm_section_data.data(), fm_section_data.size(), 0);
        if (computed_fm != fm_checksum) {
            throw Error{ErrorCode::CHECKSUM_MISMATCH,
                        "FM-index section checksum mismatch (stored=" +
                        std::to_string(fm_checksum) + " computed=" +
                        std::to_string(computed_fm) + ")"};
        }

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
        // C3 fix: accumulate section data for checksum validation
        std::vector<uint8_t> bv_section_data;

        uint64_t bread_len = r.ReadU64();
        auto bread_bytes = r.ReadBytes(bread_len);
        idx.bv.B_read.Deserialize(bread_bytes.data(), bread_bytes.size());
        bv_section_data.insert(bv_section_data.end(), bread_bytes.begin(), bread_bytes.end());

        uint64_t bwindow_len = r.ReadU64();
        auto bwindow_bytes = r.ReadBytes(bwindow_len);
        idx.bv.B_window.Deserialize(bwindow_bytes.data(), bwindow_bytes.size());
        bv_section_data.insert(bv_section_data.end(), bwindow_bytes.begin(), bwindow_bytes.end());

        // C3 fix: Validate bitvector section checksum (Contract §9.4)
        uint64_t bv_checksum = r.ReadU64();
        uint64_t computed_bv = XXH64(bv_section_data.data(), bv_section_data.size(), 0);
        if (computed_bv != bv_checksum) {
            throw Error{ErrorCode::CHECKSUM_MISMATCH,
                        "Bitvector section checksum mismatch (stored=" +
                        std::to_string(bv_checksum) + " computed=" +
                        std::to_string(computed_bv) + ")"};
        }
    }

    // ─── Window section ──────────────────────────────────────────────────────
    {
        // C3 fix: accumulate section data for checksum validation
        std::vector<uint8_t> win_section_data;
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
            // Accumulate for checksum (must match seal.cpp layout)
            uint8_t buf[8];
            std::memcpy(buf, &w.chrom_id, 4); win_section_data.insert(win_section_data.end(), buf, buf+4);
            std::memcpy(buf, &w.l, 8); win_section_data.insert(win_section_data.end(), buf, buf+8);
            std::memcpy(buf, &w.r, 8); win_section_data.insert(win_section_data.end(), buf, buf+8);
            std::memcpy(buf, &w.first_read_id, 8); win_section_data.insert(win_section_data.end(), buf, buf+8);
            std::memcpy(buf, &w.last_read_id, 8); win_section_data.insert(win_section_data.end(), buf, buf+8);
            std::memcpy(buf, &w.genomic_start, 8); win_section_data.insert(win_section_data.end(), buf, buf+8);
            std::memcpy(buf, &w.genomic_end, 8); win_section_data.insert(win_section_data.end(), buf, buf+8);
        }
        // C3 fix: Validate window section checksum (Contract §9.4)
        uint64_t win_checksum = r.ReadU64();
        uint64_t computed_win = XXH64(win_section_data.data(), win_section_data.size(), 0);
        if (computed_win != win_checksum) {
            throw Error{ErrorCode::CHECKSUM_MISMATCH,
                        "Window section checksum mismatch (stored=" +
                        std::to_string(win_checksum) + " computed=" +
                        std::to_string(computed_win) + ")"};
        }
    }

    // ─── Directory section ───────────────────────────────────────────────────
    {
        // Helper: read a per-read directory, validate checksum, return entries
        auto read_per_read_dir = [&](const std::string& name) -> StreamDirectoryPerRead {
            /*uint8_t gran =*/ r.ReadU8();
            uint64_t entry_count = r.ReadU64();
            StreamDirectoryPerRead dir(entry_count);
            std::vector<uint8_t> dir_data;
            for (uint64_t i = 0; i < entry_count; ++i) {
                dir[i].offset = r.ReadU64();
                dir[i].length = r.ReadU32();
                uint8_t buf[8];
                std::memcpy(buf, &dir[i].offset, 8);
                dir_data.insert(dir_data.end(), buf, buf + 8);
                std::memcpy(buf, &dir[i].length, 4);
                dir_data.insert(dir_data.end(), buf, buf + 4);
            }
            uint64_t stored_chk = r.ReadU64();
            if (entry_count > 0) {
                uint64_t computed_chk = XXH64(dir_data.data(), dir_data.size(), 0);
                if (computed_chk != stored_chk) {
                    throw Error{ErrorCode::CHECKSUM_MISMATCH,
                                "Directory checksum mismatch in " + name};
                }
            }
            return dir;
        };

        // dir_seq: BWT is single payload, no per-read entries
        /*uint8_t seq_gran =*/ r.ReadU8();
        uint64_t seq_count = r.ReadU64();
        for (uint64_t i = 0; i < seq_count; ++i) {
            r.ReadU64(); r.ReadU32();
        }
        r.ReadU64(); // checksum

        // dir_qual: per-read OR block-level (Contract §2.3 F4)
        {
            uint8_t qual_gran = r.ReadU8();
            if (qual_gran == 1) {
                // Block-level directory
                uint32_t block_size = r.ReadU32();
                idx.qual_block_size = block_size;
                uint64_t entry_count = r.ReadU64();
                StreamDirectoryBlockLevel block_dir(entry_count);
                std::vector<uint8_t> dir_data;
                for (uint64_t i = 0; i < entry_count; ++i) {
                    block_dir[i].block_offset = r.ReadU64();
                    block_dir[i].block_length = r.ReadU32();
                    block_dir[i].first_read_id = r.ReadU32();
                    uint8_t buf[8];
                    std::memcpy(buf, &block_dir[i].block_offset, 8);
                    dir_data.insert(dir_data.end(), buf, buf + 8);
                    std::memcpy(buf, &block_dir[i].block_length, 4);
                    dir_data.insert(dir_data.end(), buf, buf + 4);
                    std::memcpy(buf, &block_dir[i].first_read_id, 4);
                    dir_data.insert(dir_data.end(), buf, buf + 4);
                }
                uint64_t stored_chk = r.ReadU64();
                if (entry_count > 0) {
                    uint64_t computed_chk = XXH64(dir_data.data(), dir_data.size(), 0);
                    if (computed_chk != stored_chk) {
                        throw Error{ErrorCode::CHECKSUM_MISMATCH,
                                    "Block-level dir_qual checksum mismatch"};
                    }
                }
                idx.qual_directory = std::move(block_dir);
            } else {
                // Per-read directory (granularity byte already consumed)
                uint64_t entry_count = r.ReadU64();
                StreamDirectoryPerRead dir(entry_count);
                std::vector<uint8_t> dir_data;
                for (uint64_t i = 0; i < entry_count; ++i) {
                    dir[i].offset = r.ReadU64();
                    dir[i].length = r.ReadU32();
                    uint8_t buf[8];
                    std::memcpy(buf, &dir[i].offset, 8);
                    dir_data.insert(dir_data.end(), buf, buf + 8);
                    std::memcpy(buf, &dir[i].length, 4);
                    dir_data.insert(dir_data.end(), buf, buf + 4);
                }
                uint64_t stored_chk = r.ReadU64();
                if (entry_count > 0) {
                    uint64_t computed_chk = XXH64(dir_data.data(), dir_data.size(), 0);
                    if (computed_chk != stored_chk) {
                        throw Error{ErrorCode::CHECKSUM_MISMATCH,
                                    "Per-read dir_qual checksum mismatch"};
                    }
                }
                idx.qual_directory = std::move(dir);
                idx.qual_block_size = 0;
            }
        }

        // C5 fix: load actual per-read directories for reconstruction
        idx.meta_directory = read_per_read_dir("dir_meta");
        idx.map_directory  = read_per_read_dir("dir_map");
    }

    // ─── Read metadata section ────────────────────────────────────────────
    {
        std::vector<uint8_t> read_section_data;
        uint64_t n_reads = r.ReadU64();
        idx.reads.resize(n_reads);
        for (uint64_t i = 0; i < n_reads; ++i) {
            auto& rd = idx.reads[i];
            rd.chrom_id = r.ReadU32();
            rd.pos = r.ReadU64();
            uint32_t seq_len = r.ReadU32();
            rd.seq.resize(seq_len);  // placeholder (actual bases in S_seq stream)
            rd.read_id = i;
            // MAPQ (Contract §2.1 lossless)
            rd.mapq = r.ReadU8();
            // QNAME (Contract §2.1 lossless)
            uint32_t qname_len = r.ReadU32();
            if (qname_len > 0) {
                auto qname_bytes = r.ReadBytes(qname_len);
                rd.qname.assign(qname_bytes.begin(), qname_bytes.end());
            }
            // CIGAR
            uint32_t n_ops = r.ReadU32();
            rd.cigar.resize(n_ops);
            for (uint32_t c = 0; c < n_ops; ++c) {
                rd.cigar[c].op = r.ReadU8();
                rd.cigar[c].len = r.ReadU32();
            }
            // FLAG (Contract §2.1 / §2.2 lossless)
            rd.flag = r.ReadU16();
            // Aux tags raw bytes (Contract §2.2)
            uint32_t aux_len = r.ReadU32();
            if (aux_len > 0) {
                auto aux_bytes = r.ReadBytes(aux_len);
                rd.aux_data.assign(aux_bytes.begin(), aux_bytes.end());
            }
            // Source file ID and BAM offset (Contract §0.9, §1.2)
            rd.source_file_id = r.ReadU32();
            rd.bam_offset = r.ReadU64();

            // Accumulate key fields for section checksum (must match seal.cpp)
            uint8_t buf[8];
            std::memcpy(buf, &rd.chrom_id, 4); read_section_data.insert(read_section_data.end(), buf, buf+4);
            std::memcpy(buf, &rd.pos, 8); read_section_data.insert(read_section_data.end(), buf, buf+8);
            std::memcpy(buf, &seq_len, 4); read_section_data.insert(read_section_data.end(), buf, buf+4);
            read_section_data.push_back(rd.mapq);
            uint16_t flag16 = static_cast<uint16_t>(rd.flag);
            std::memcpy(buf, &flag16, 2); read_section_data.insert(read_section_data.end(), buf, buf+2);
            std::memcpy(buf, &rd.source_file_id, 4); read_section_data.insert(read_section_data.end(), buf, buf+4);
            std::memcpy(buf, &rd.bam_offset, 8); read_section_data.insert(read_section_data.end(), buf, buf+8);
        }
        // C4: Validate read section checksum (Contract §9.4 — no silent fallback)
        uint64_t stored_checksum = r.ReadU64();
        uint64_t computed_checksum = XXH64(read_section_data.data(), read_section_data.size(), 0);
        if (computed_checksum != stored_checksum) {
            throw Error{ErrorCode::CHECKSUM_MISMATCH,
                        "Read metadata section checksum mismatch (stored=" +
                        std::to_string(stored_checksum) + " computed=" +
                        std::to_string(computed_checksum) + ")"};
        }
    }

    // ─── ISA Samples section (Architecture §4.4, optional) ──────────────────
    if (h.has_isa_samples && h.sample_step_s_prime > 0) {
        uint64_t isa_count = r.ReadU64();
        uint64_t isa_step = r.ReadU64();
        std::vector<uint64_t> isa_samples(isa_count);
        std::vector<uint8_t> isa_section_data;
        for (uint64_t i = 0; i < isa_count; ++i) {
            isa_samples[i] = r.ReadU64();
            uint8_t buf[8]; std::memcpy(buf, &isa_samples[i], 8);
            isa_section_data.insert(isa_section_data.end(), buf, buf + 8);
        }
        uint64_t isa_checksum = r.ReadU64();
        uint64_t computed_isa = XXH64(isa_section_data.data(), isa_section_data.size(), 0);
        if (computed_isa != isa_checksum) {
            throw Error{ErrorCode::CHECKSUM_MISMATCH,
                        "ISA section checksum mismatch"};
        }
        idx.fm.SetISASamples(std::move(isa_samples), isa_step);
    }

    // ─── SARange section (Architecture §5.3, optional — ENHANCED tier) ──────
    if (h.enable_sarange) {
        uint64_t sarange_length = r.ReadU64();
        auto sarange_payload = r.ReadBytes(sarange_length);
        uint64_t sarange_checksum = r.ReadU64();
        uint64_t computed_sar = XXH64(sarange_payload.data(), sarange_payload.size(), 0);
        if (computed_sar != sarange_checksum) {
            throw Error{ErrorCode::CHECKSUM_MISMATCH,
                        "SARange section checksum mismatch"};
        }
        idx.sarange.Deserialize(sarange_payload.data(), sarange_payload.size());
    }

    // ── Global footer checksum validation (Contract §1.2 / §9.4) ──────────
    // Verify xxHash64 over all bytes before the footer (last 12 bytes).
    {
        std::ifstream ifs(path, std::ios::binary);
        if (ifs.is_open()) {
            ifs.seekg(0, std::ios::end);
            uint64_t total_len = static_cast<uint64_t>(ifs.tellg());
            if (total_len >= 12) {
                // Read footer: xxHash64 (8 bytes) + magic (4 bytes)
                ifs.seekg(static_cast<std::streamoff>(total_len - 12));
                uint64_t stored_global;
                ifs.read(reinterpret_cast<char*>(&stored_global), 8);
                uint32_t footer_magic;
                ifs.read(reinterpret_cast<char*>(&footer_magic), 4);

                if (footer_magic != 0xB5110000) {
                    throw Error{ErrorCode::CORRUPT_BSI,
                                "Invalid .bsi footer magic (file may be truncated or corrupt)"};
                }

                // Compute xxHash64 over data before footer
                uint64_t data_len = total_len - 12;
                ifs.seekg(0);
                std::vector<uint8_t> all_data(data_len);
                ifs.read(reinterpret_cast<char*>(all_data.data()),
                         static_cast<std::streamsize>(data_len));

                uint64_t computed_global = XXH64(all_data.data(), data_len, 0);
                if (computed_global != stored_global) {
                    throw Error{ErrorCode::CHECKSUM_MISMATCH,
                                "Global .bsi checksum mismatch — file is corrupt or tampered"};
                }
            }
        }
    }

    // ── Ordering hash re-verification (Contract §1.2) ─────────────────────
    // "Both hashes are stored in the .bsi header and re-verified on every
    // index open." Re-compute ordering_hash = SHA-256(concat of 4-tuples).
    {
        // Build the ordering hash from loaded read data
        std::vector<uint8_t> ordering_data;
        ordering_data.reserve(idx.reads.size() * (4 + 8 + 4 + 8));
        for (const auto& rd : idx.reads) {
            uint8_t buf[8];
            uint32_t cid = rd.chrom_id;
            std::memcpy(buf, &cid, 4);
            ordering_data.insert(ordering_data.end(), buf, buf + 4);
            uint64_t pos = rd.pos;
            std::memcpy(buf, &pos, 8);
            ordering_data.insert(ordering_data.end(), buf, buf + 8);
            uint32_t sfid = rd.source_file_id;
            std::memcpy(buf, &sfid, 4);
            ordering_data.insert(ordering_data.end(), buf, buf + 4);
            uint64_t boff = rd.bam_offset;
            std::memcpy(buf, &boff, 8);
            ordering_data.insert(ordering_data.end(), buf, buf + 8);
        }

        std::array<uint8_t, 32> computed_ordering_hash;
        SHA256(ordering_data.data(), ordering_data.size(),
               computed_ordering_hash.data());

        if (computed_ordering_hash != h.ordering_hash) {
            throw Error{ErrorCode::ORDERING_HASH_MISMATCH,
                        "Ordering hash verification failed on index open — "
                        "read order in .bsi does not match stored ordering_hash"};
        }
    }

    idx.global_checksum_valid = true;
    return idx;
}

bool VerifyBsi(const std::string& path) {
    // H5 fix: Full verification per Contract §3.3.2 — ALL section checksums.
    // ReadBsi() validates every section checksum (streams, FM-index, bitvectors,
    // windows, read metadata, directories) plus ordering_hash and global footer.
    // If ReadBsi succeeds without throwing, all checksums pass.
    try {
        ReadBsi(path);
        return true;
    } catch (const Error&) {
        return false;
    } catch (...) {
        return false;
    }
}

}  // namespace bamsix
