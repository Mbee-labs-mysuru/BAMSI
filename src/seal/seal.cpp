#include "seal.hpp"

#define XXH_INLINE_ALL
#include <xxhash.h>

#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

namespace bamsi {

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
        uint32_t footer_magic = 0xB5110000;  // BAMSI footer marker
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
    w.Write(h.bamsi_version, 16);
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
        auto bwt_bytes = input.fm.SerializeBWT();
        uint64_t bwt_length = bwt_bytes.size();
        w.WriteU64(bwt_length);
        w.WriteBytes(bwt_bytes);

        // C array: 7 × uint64
        for (auto val : input.fm.CArray()) {
            w.WriteU64(val);
        }

        // Occ metadata
        auto occ_bytes = input.fm.Occ().Serialize();
        uint64_t occ_len = occ_bytes.size();
        w.WriteU64(occ_len);
        w.WriteBytes(occ_bytes);

        // SA samples
        const auto& sa_samples = input.fm.SASamples();
        uint64_t sa_count = sa_samples.size();
        w.WriteU64(sa_count);
        for (auto val : sa_samples) {
            w.WriteU64(val);
        }

        // Section checksum placeholder (simplified)
        w.WriteU64(0);  // TODO: proper section checksum
    }

    // ─── Bitvector section ───────────────────────────────────────────────────
    {
        auto bread = input.bv.B_read.Serialize();
        uint64_t bread_len = bread.size();
        w.WriteU64(bread_len);
        w.WriteBytes(bread);

        auto bwindow = input.bv.B_window.Serialize();
        uint64_t bwindow_len = bwindow.size();
        w.WriteU64(bwindow_len);
        w.WriteBytes(bwindow);

        w.WriteU64(0);  // section checksum placeholder
    }

    // ─── Window section ──────────────────────────────────────────────────────
    {
        for (const auto& win : input.windows) {
            w.WriteU32(win.chrom_id);
            w.WriteU64(win.l);
            w.WriteU64(win.r);
            w.WriteU64(win.first_read_id);
            w.WriteU64(win.last_read_id);
            w.WriteU64(win.genomic_start);
            w.WriteU64(win.genomic_end);
        }
        w.WriteU64(0);  // section checksum placeholder
    }

    // ─── Directory section ───────────────────────────────────────────────────
    {
        // dir_seq (per-read for V1)
        auto write_per_read_dir = [&](const StreamDirectoryPerRead& dir) {
            w.WriteU8(0);  // granularity = PER_READ
            uint64_t entry_count = dir.size();
            w.WriteU64(entry_count);
            for (const auto& e : dir) {
                w.WriteU64(e.offset);
                w.WriteU32(e.length);
            }
            w.WriteU64(0);  // section checksum placeholder
        };

        // For V1, dir_seq and dir_qual are empty per-read directories
        // (BWT is stored as single payload, not per-read).
        // We write empty directories for these.
        w.WriteU8(0);  // dir_seq granularity
        w.WriteU64(0); // dir_seq entry count = 0 (single payload)
        w.WriteU64(0); // checksum

        w.WriteU8(0);  // dir_qual granularity
        w.WriteU64(0); // dir_qual entry count = 0 placeholder
        w.WriteU64(0); // checksum

        // dir_meta: mandatory per-read
        write_per_read_dir(input.meta.directory);

        // dir_map: mandatory per-read
        write_per_read_dir(input.map.directory);
    }

    // ─── Read metadata section (for Locate at query time) ─────────────────
    {
        uint64_t n_reads = input.reads.size();
        w.WriteU64(n_reads);
        for (const auto& rd : input.reads) {
            w.WriteU32(rd.chrom_id);
            w.WriteU64(rd.pos);
            uint32_t seq_len = static_cast<uint32_t>(rd.seq.size());
            w.WriteU32(seq_len);
            // CIGAR
            uint32_t n_ops = static_cast<uint32_t>(rd.cigar.size());
            w.WriteU32(n_ops);
            for (const auto& cop : rd.cigar) {
                w.WriteU8(cop.op);
                w.WriteU32(cop.len);
            }
        }
        w.WriteU64(0);  // section checksum placeholder
    }

    // ─── ISA Samples section (Architecture §4.4, optional) ──────────────────
    if (!input.isa_samples.empty() && input.isa_step > 0) {
        uint64_t isa_count = input.isa_samples.size();
        w.WriteU64(isa_count);
        w.WriteU64(input.isa_step);
        for (auto val : input.isa_samples) {
            w.WriteU64(val);
        }
        w.WriteU64(0);  // section checksum placeholder
    }

    // ─── SARange section (Architecture §5.3, optional — ENHANCED tier) ──────
    if (input.sarange.IsBuilt()) {
        auto sarange_payload = input.sarange.Serialize();
        uint64_t sarange_length = sarange_payload.size();
        w.WriteU64(sarange_length);
        w.WriteBytes(sarange_payload);
        w.WriteU64(0);  // section checksum placeholder
    }

    // ─── Footer ──────────────────────────────────────────────────────────────
    w.Finalize();
}

}  // namespace bamsi
