#include "ingest.hpp"

#include <htslib/sam.h>
#include <htslib/hts.h>

#include <algorithm>
#include <cstring>
#include <iostream>
#include <set>
#include <openssl/sha.h>

namespace bamsi {

namespace {

/// Convert a 4-bit BAM base code to our 0..4 code.
/// BAM uses: =0, A1, C2, M3, G4, R5, S6, V7, T8, W9, YA, HB, KC, DD, NF
uint8_t BamBaseToCode(uint8_t bam4bit) {
    switch (bam4bit) {
        case 1:  return CODE_A;  // A
        case 2:  return CODE_C;  // C
        case 4:  return CODE_G;  // G
        case 8:  return CODE_T;  // T
        default: return CODE_N;  // N or ambiguous
    }
}

/// Compute SHA-256 of a byte buffer.
std::array<uint8_t, 32> Sha256(const uint8_t* data, size_t len) {
    std::array<uint8_t, 32> digest = {};
    SHA256(data, len, digest.data());
    return digest;
}

/// Compute SHA-256 of a vector.
std::array<uint8_t, 32> Sha256(const std::vector<uint8_t>& data) {
    return Sha256(data.data(), data.size());
}

}  // namespace

IngestResult IngestBams(const std::vector<std::string>& bam_paths) {
    IngestResult result;
    result.source_file_count = static_cast<uint32_t>(bam_paths.size());

    // Collect all distinct chromosome names across all files
    std::set<std::string> all_chroms;
    // Per-file header bytes for manifest hash
    std::vector<std::vector<uint8_t>> header_bytes_per_file;

    // First pass: collect chromosome names and header bytes
    for (uint32_t file_id = 0; file_id < bam_paths.size(); ++file_id) {
        htsFile* fp = hts_open(bam_paths[file_id].c_str(), "r");
        if (!fp) {
            throw Error{ErrorCode::INVALID_BAM_INPUT,
                        "Cannot open BAM file: " + bam_paths[file_id]};
        }

        sam_hdr_t* hdr = sam_hdr_read(fp);
        if (!hdr) {
            hts_close(fp);
            throw Error{ErrorCode::INVALID_BAM_INPUT,
                        "Cannot read BAM header: " + bam_paths[file_id]};
        }

        // Collect reference names
        for (int i = 0; i < sam_hdr_nref(hdr); ++i) {
            all_chroms.insert(sam_hdr_tid2name(hdr, i));
        }

        // Store raw header text bytes for manifest hash
        const char* hdr_text = sam_hdr_str(hdr);
        int hdr_len = sam_hdr_length(hdr);
        std::vector<uint8_t> hdr_bytes(hdr_text, hdr_text + hdr_len);
        header_bytes_per_file.push_back(std::move(hdr_bytes));

        sam_hdr_destroy(hdr);
        hts_close(fp);
    }

    // Build frozen chrom_id mapping: lexicographic sort of all distinct names
    result.chrom_names.assign(all_chroms.begin(), all_chroms.end());
    std::sort(result.chrom_names.begin(), result.chrom_names.end());
    for (uint32_t i = 0; i < result.chrom_names.size(); ++i) {
        result.chrom_to_id[result.chrom_names[i]] = i;
    }

    // Compute source_manifest_hash per §0.10:
    // SHA-256(concat for f=0..F-1: uint32_le(len(filename_f)) || utf8(filename_f) || SHA-256(header_f))
    {
        SHA256_CTX ctx;
        SHA256_Init(&ctx);
        for (uint32_t f = 0; f < bam_paths.size(); ++f) {
            const std::string& fname = bam_paths[f];
            uint32_t name_len = static_cast<uint32_t>(fname.size());
            // uint32_le
            SHA256_Update(&ctx, &name_len, sizeof(name_len));
            // utf8 bytes
            SHA256_Update(&ctx, fname.data(), fname.size());
            // SHA-256 of header bytes (32-byte digest)
            auto hdr_hash = Sha256(header_bytes_per_file[f]);
            SHA256_Update(&ctx, hdr_hash.data(), 32);
        }
        SHA256_Final(result.source_manifest_hash.data(), &ctx);
    }

    // Second pass: read all records, apply inclusion rule, extract RawReads
    for (uint32_t file_id = 0; file_id < bam_paths.size(); ++file_id) {
        htsFile* fp = hts_open(bam_paths[file_id].c_str(), "r");
        sam_hdr_t* hdr = sam_hdr_read(fp);
        bam1_t* rec = bam_init1();

        uint64_t record_index = 0;  // bam_offset: sequential record index in file

        while (sam_read1(fp, hdr, rec) >= 0) {
            uint16_t flag = rec->core.flag;

            // Inclusion rule (§0.1):
            // Mapped, primary, non-supplementary, POS >= 1
            bool is_mapped         = !(flag & BAM_FUNMAP);         // FLAG & 0x4 == 0
            bool is_primary        = !(flag & BAM_FSECONDARY);     // FLAG & 0x100 == 0
            bool is_not_supp       = !(flag & BAM_FSUPPLEMENTARY); // FLAG & 0x800 == 0
            // BAM POS is 0-based; SAM POS is 1-based. rec->core.pos is 0-based.
            bool has_valid_pos     = (rec->core.pos >= 0);         // POS >= 1 in 1-based

            if (is_mapped && is_primary && is_not_supp && has_valid_pos) {
                RawRead read;

                // Extract sequence
                int32_t seq_len = rec->core.l_qseq;
                read.seq.resize(seq_len);
                uint8_t* bam_seq = bam_get_seq(rec);
                for (int32_t i = 0; i < seq_len; ++i) {
                    read.seq[i] = BamBaseToCode(bam_seqi(bam_seq, i));
                }

                // Chromosome name
                int tid = rec->core.tid;
                if (tid >= 0) {
                    read.chrom = sam_hdr_tid2name(hdr, tid);
                }

                // FLAG
                read.flag = flag;

                // POS: convert BAM 0-based to SAM 1-based (exactly once)
                read.pos = static_cast<uint64_t>(rec->core.pos) + 1;

                // CIGAR
                uint32_t* cigar_raw = bam_get_cigar(rec);
                uint32_t n_cigar = rec->core.n_cigar;
                read.cigar.resize(n_cigar);
                for (uint32_t c = 0; c < n_cigar; ++c) {
                    read.cigar[c].op  = static_cast<uint8_t>(bam_cigar_op(cigar_raw[c]));
                    read.cigar[c].len = bam_cigar_oplen(cigar_raw[c]);
                }

                // Quality scores
                uint8_t* qual_ptr = bam_get_qual(rec);
                read.qual.resize(seq_len);
                for (int32_t i = 0; i < seq_len; ++i) {
                    read.qual[i] = (qual_ptr[i] == 0xFF) ? 0 : qual_ptr[i];
                }

                // Verify |qual| == |seq|
                if (read.qual.size() != read.seq.size()) {
                    throw Error{ErrorCode::INVALID_BAM_INPUT,
                                "Quality length mismatch for record"};
                }

                read.source_file_id = file_id;
                read.bam_offset = record_index;

                result.reads.push_back(std::move(read));
            }
            ++record_index;
        }

        bam_destroy1(rec);
        sam_hdr_destroy(hdr);
        hts_close(fp);
    }

    return result;
}

}  // namespace bamsi
