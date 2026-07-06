#include "seal.hpp"

#define XXH_INLINE_ALL
#include <xxhash.h>

#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

namespace bamsix {

namespace {

/// Write raw bytes to output and accumulate xxHash.
class BsiWriter {
public:
    explicit BsiWriter(const std::string& path)
        : ofs_(path, std::ios::binary) {
        if (!ofs_.is_open()) {
            throw Error{ErrorCode::BUILD_VALIDATION_FAILED,
                        "Cannot open output file: " + path};
        }
        global_hash_ = XXH64_createState();
        XXH64_reset(global_hash_, 0);
    }

    ~BsiWriter() {
        XXH64_freeState(global_hash_);
    }

    void Write(const void* data, size_t len) {
        ofs_.write(static_cast<const char*>(data), len);
        XXH64_update(global_hash_, data, len);
        bytes_written_ += len;
    }

    void WriteU8(uint8_t val) { Write(&val, 1); }
    void WriteU16(uint16_t val) { Write(&val, 2); }
    void WriteU32(uint32_t val) { Write(&val, 4); }
    void WriteU64(uint64_t val) { Write(&val, 8); }
    void WriteBytes(const std::vector<uint8_t>& data) {
        Write(data.data(), data.size());
    }
    void WriteBytes(const uint8_t* data, size_t len) {
        Write(data, len);
    }

    /// Compute section checksum (xxHash64 of given bytes).
    static uint64_t SectionChecksum(const uint8_t* data, size_t len) {
        return XXH64(data, len, 0);
    }

    /// Write a stream section: payload_length, codec_id, payload, section_checksum.
    void WriteStreamSection(const std::vector<uint8_t>& payload, uint8_t codec_id) {
        uint64_t payload_length = payload.size();
        WriteU64(payload_length);
        WriteU8(codec_id);
        WriteBytes(payload);

        // Section checksum: xxHash64(codec_id || payload)
        std::vector<uint8_t> check_data;
        check_data.push_back(codec_id);
        check_data.insert(check_data.end(), payload.begin(), payload.end());
        uint64_t checksum = SectionChecksum(check_data.data(), check_data.size());
        WriteU64(checksum);
    }

    /// Finalize: write global checksum and footer magic.
    void Finalize() {
        uint64_t global_checksum = XXH64_digest(global_hash_);
        ofs_.write(reinterpret_cast<const char*>(&global_checksum), 8);
        uint32_t footer_magic = 0xB5110000;  // BAMSIX footer marker
        ofs_.write(reinterpret_cast<const char*>(&footer_magic), 4);
        ofs_.close();
    }

    uint64_t BytesWritten() const { return bytes_written_; }

private:
    std::ofstream ofs_;
    XXH64_state_t* global_hash_;
    uint64_t bytes_written_ = 0;
};

}  // namespace

void WriteBsi(const std::string& path, const SealInput& input) {
    BsiWriter w(path);
    const auto& h = input.header;

    // ─── Header ──────────────────────────────────────────────────────────────
    w.Write(h.magic, 4);
    w.WriteU16(h.version);
    w.Write(h.bamsix_version, 16);
    w.WriteU8(h.host_os_id);
    w.WriteU8(h.cpu_arch_id);
    w.WriteU64(h.build_timestamp_utc);
    w.WriteU8(h.is_lossless);
    w.WriteU32(h.source_file_count);
    w.WriteBytes(h.source_manifest_hash.data(), 32);
    w.WriteBytes(h.ordering_hash.data(), 32);
    w.WriteU64(h.S_length);
    w.WriteU64(h.N_reads);
    w.WriteU32(h.N_windows);
    w.WriteU32(h.sample_step_s);
    w.WriteU8(h.has_isa_samples);
    w.WriteU32(h.sample_step_s_prime);
    w.WriteU8(h.enable_sarange);
    w.WriteU8(h.sarange_variant);
    w.WriteU8(h.shared_bwt);
    w.WriteU8(h.enable_bidirectional);
    w.WriteU8(h.recommended_seed_length);
    w.WriteU64(h.window_size_T);
    w.WriteU8(h.entropy_order_k);
    w.WriteU8(h.qual_codec_id);
    w.WriteU8(h.qual_lossy_bins);
    w.WriteU8(h.meta_codec_id);
    w.WriteU8(h.map_codec_id);
    w.WriteU8(h.strand_mode);
    w.WriteU64(h.sentinel_row);
    w.WriteU32(h.chrom_count);

    // Chromosome name table
    for (const auto& ce : h.chrom_name_table) {
        w.WriteU32(ce.chrom_id);
        uint32_t name_len = static_cast<uint32_t>(ce.name.size());
        w.WriteU32(name_len);
        w.Write(ce.name.data(), ce.name.size());
    }

    w.WriteU32(h.seq_block_size);
    w.WriteU32(h.qual_block_size);
    w.WriteU8(h.allow_parallel_sa);
    w.WriteU8(h.reference_based_encoding);
    w.WriteBytes(h.reference_sha256.data(), 32);
    w.WriteU32(h.flags);

    // ─── Stream sections ─────────────────────────────────────────────────────
    w.WriteStreamSection(input.seq.payload, input.seq.codec_id);
    w.WriteStreamSection(input.qual.payload, static_cast<uint8_t>(input.qual.codec_id));
    w.WriteStreamSection(input.meta.payload, static_cast<uint8_t>(input.meta.codec_id));
    w.WriteStreamSection(input.map.payload, static_cast<uint8_t>(input.map.codec_id));

    // ─── FM-Index section ────────────────────────────────────────────────────
    {
        // Accumulate section bytes for checksum computation
        std::vector<uint8_t> fm_section_data;

        auto bwt_bytes = input.fm.SerializeBWT();
        uint64_t bwt_length = bwt_bytes.size();
        w.WriteU64(bwt_length);
        w.WriteBytes(bwt_bytes);
        fm_section_data.insert(fm_section_data.end(), bwt_bytes.begin(), bwt_bytes.end());

        // C array: 7 × uint64
        for (auto val : input.fm.CArray()) {
            w.WriteU64(val);
            uint8_t buf[8]; std::memcpy(buf, &val, 8);
            fm_section_data.insert(fm_section_data.end(), buf, buf + 8);
        }

        // Occ metadata
        auto occ_bytes = input.fm.Occ().Serialize();
        uint64_t occ_len = occ_bytes.size();
        w.WriteU64(occ_len);
        w.WriteBytes(occ_bytes);
        fm_section_data.insert(fm_section_data.end(), occ_bytes.begin(), occ_bytes.end());

        // SA samples
        const auto& sa_samples = input.fm.SASamples();
        uint64_t sa_count = sa_samples.size();
        w.WriteU64(sa_count);
        for (auto val : sa_samples) {
            w.WriteU64(val);
            uint8_t buf[8]; std::memcpy(buf, &val, 8);
            fm_section_data.insert(fm_section_data.end(), buf, buf + 8);
        }

        // Section checksum (Contract §9.4)
        uint64_t fm_checksum = XXH64(fm_section_data.data(), fm_section_data.size(), 0);
        w.WriteU64(fm_checksum);
    }

    // ─── Bitvector section ───────────────────────────────────────────────────
    {
        std::vector<uint8_t> bv_section_data;

        auto bread = input.bv.B_read.Serialize();
        uint64_t bread_len = bread.size();
        w.WriteU64(bread_len);
        w.WriteBytes(bread);
        bv_section_data.insert(bv_section_data.end(), bread.begin(), bread.end());

        auto bwindow = input.bv.B_window.Serialize();
        uint64_t bwindow_len = bwindow.size();
        w.WriteU64(bwindow_len);
        w.WriteBytes(bwindow);
        bv_section_data.insert(bv_section_data.end(), bwindow.begin(), bwindow.end());

        uint64_t bv_checksum = XXH64(bv_section_data.data(), bv_section_data.size(), 0);
        w.WriteU64(bv_checksum);
    }

    // ─── Window section ──────────────────────────────────────────────────────
    {
        std::vector<uint8_t> win_section_data;
        for (const auto& win : input.windows) {
            w.WriteU32(win.chrom_id);
            w.WriteU64(win.l);
            w.WriteU64(win.r);
            w.WriteU64(win.first_read_id);
            w.WriteU64(win.last_read_id);
            w.WriteU64(win.genomic_start);
            w.WriteU64(win.genomic_end);
            // Accumulate for checksum
            uint8_t buf[8];
            std::memcpy(buf, &win.chrom_id, 4); win_section_data.insert(win_section_data.end(), buf, buf+4);
            std::memcpy(buf, &win.l, 8); win_section_data.insert(win_section_data.end(), buf, buf+8);
            std::memcpy(buf, &win.r, 8); win_section_data.insert(win_section_data.end(), buf, buf+8);
            std::memcpy(buf, &win.first_read_id, 8); win_section_data.insert(win_section_data.end(), buf, buf+8);
            std::memcpy(buf, &win.last_read_id, 8); win_section_data.insert(win_section_data.end(), buf, buf+8);
            std::memcpy(buf, &win.genomic_start, 8); win_section_data.insert(win_section_data.end(), buf, buf+8);
            std::memcpy(buf, &win.genomic_end, 8); win_section_data.insert(win_section_data.end(), buf, buf+8);
        }
        uint64_t win_checksum = XXH64(win_section_data.data(), win_section_data.size(), 0);
        w.WriteU64(win_checksum);
    }

    // ─── Directory section ───────────────────────────────────────────────────
    {
        // H3 fix: compute real checksums for every directory section
        auto write_per_read_dir = [&](const StreamDirectoryPerRead& dir) {
            std::vector<uint8_t> dir_data;
            w.WriteU8(0);  // granularity = PER_READ
            uint64_t entry_count = dir.size();
            w.WriteU64(entry_count);
            for (const auto& e : dir) {
                w.WriteU64(e.offset);
                w.WriteU32(e.length);
                uint8_t buf[8];
                std::memcpy(buf, &e.offset, 8);
                dir_data.insert(dir_data.end(), buf, buf + 8);
                std::memcpy(buf, &e.length, 4);
                dir_data.insert(dir_data.end(), buf, buf + 4);
            }
            uint64_t dir_checksum = BsiWriter::SectionChecksum(dir_data.data(), dir_data.size());
            w.WriteU64(dir_checksum);
        };

        // Contract §2.3 F4: block-level directory writer
        auto write_block_level_dir = [&](const StreamDirectoryBlockLevel& dir, uint32_t block_size) {
            std::vector<uint8_t> dir_data;
            w.WriteU8(1);  // granularity = BLOCK_LEVEL
            w.WriteU32(block_size);  // B_dir stored so decoders can locate blocks
            uint64_t entry_count = dir.size();
            w.WriteU64(entry_count);
            for (const auto& e : dir) {
                w.WriteU64(e.block_offset);
                w.WriteU32(e.block_length);
                w.WriteU32(e.first_read_id);
                uint8_t buf[8];
                std::memcpy(buf, &e.block_offset, 8);
                dir_data.insert(dir_data.end(), buf, buf + 8);
                std::memcpy(buf, &e.block_length, 4);
                dir_data.insert(dir_data.end(), buf, buf + 4);
                std::memcpy(buf, &e.first_read_id, 4);
                dir_data.insert(dir_data.end(), buf, buf + 4);
            }
            uint64_t dir_checksum = BsiWriter::SectionChecksum(dir_data.data(), dir_data.size());
            w.WriteU64(dir_checksum);
        };

        // dir_seq: BWT is stored as single payload, no per-read directory
        w.WriteU8(0);  // dir_seq granularity
        w.WriteU64(0); // dir_seq entry count = 0 (single payload)
        w.WriteU64(0); // checksum over empty data = 0

        // dir_qual: per-read OR block-level depending on qual_block_size (F4)
        if (std::holds_alternative<StreamDirectoryBlockLevel>(input.qual.directory)) {
            write_block_level_dir(std::get<StreamDirectoryBlockLevel>(input.qual.directory),
                                  input.qual.block_size);
        } else {
            write_per_read_dir(std::get<StreamDirectoryPerRead>(input.qual.directory));
        }

        // dir_meta: mandatory per-read
        write_per_read_dir(input.meta.directory);

        // dir_map: mandatory per-read
        write_per_read_dir(input.map.directory);
    }

    // ─── Read metadata section (for Locate at query time) ─────────────────
    {
        std::vector<uint8_t> read_section_data;
        uint64_t n_reads = input.reads.size();
        w.WriteU64(n_reads);
        for (const auto& rd : input.reads) {
            w.WriteU32(rd.chrom_id);
            w.WriteU64(rd.pos);
            uint32_t seq_len = static_cast<uint32_t>(rd.seq.size());
            w.WriteU32(seq_len);
            // MAPQ (Contract §2.1 lossless)
            w.WriteU8(rd.mapq);
            // QNAME length + bytes (Contract §2.1 lossless)
            uint32_t qname_len = static_cast<uint32_t>(rd.qname.size());
            w.WriteU32(qname_len);
            if (qname_len > 0) {
                w.Write(rd.qname.data(), qname_len);
            }
            // CIGAR
            uint32_t n_ops = static_cast<uint32_t>(rd.cigar.size());
            w.WriteU32(n_ops);
            for (const auto& cop : rd.cigar) {
                w.WriteU8(cop.op);
                w.WriteU32(cop.len);
            }
            // FLAG (Contract §2.1 / §2.2 lossless — S_meta includes FLAG)
            w.WriteU16(static_cast<uint16_t>(rd.flag));
            // Aux tags raw bytes (Contract §2.2 — S_meta includes optional BAM tags)
            uint32_t aux_len = static_cast<uint32_t>(rd.aux_data.size());
            w.WriteU32(aux_len);
            if (aux_len > 0) {
                w.Write(rd.aux_data.data(), aux_len);
            }
            // Source file ID and BAM offset (Contract §0.9, §1.2 — for ordering hash re-verification)
            w.WriteU32(rd.source_file_id);
            w.WriteU64(rd.bam_offset);
            // Accumulate ALL key fields for section checksum
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
        uint64_t read_checksum = XXH64(read_section_data.data(), read_section_data.size(), 0);
        w.WriteU64(read_checksum);
    }

    // ─── ISA Samples section (Architecture §4.4, optional) ──────────────────
    if (!input.isa_samples.empty() && input.isa_step > 0) {
        // H2 fix: compute real checksum over ISA section data
        std::vector<uint8_t> isa_section_data;
        uint64_t isa_count = input.isa_samples.size();
        w.WriteU64(isa_count);
        w.WriteU64(input.isa_step);
        for (auto val : input.isa_samples) {
            w.WriteU64(val);
            uint8_t buf[8]; std::memcpy(buf, &val, 8);
            isa_section_data.insert(isa_section_data.end(), buf, buf + 8);
        }
        uint64_t isa_checksum = XXH64(isa_section_data.data(), isa_section_data.size(), 0);
        w.WriteU64(isa_checksum);
    }

    // ─── SARange section (Architecture §5.3, optional — ENHANCED tier) ──────
    if (input.sarange.IsBuilt()) {
        auto sarange_payload = input.sarange.Serialize();
        uint64_t sarange_length = sarange_payload.size();
        w.WriteU64(sarange_length);
        w.WriteBytes(sarange_payload);
        // H2 fix: compute real checksum over SARange section data
        uint64_t sarange_checksum = XXH64(sarange_payload.data(), sarange_payload.size(), 0);
        w.WriteU64(sarange_checksum);
    }

    // ─── Footer ──────────────────────────────────────────────────────────────
    w.Finalize();
}

}  // namespace bamsix
