#include "build.hpp"

#include "../ingest/ingest.hpp"
#include "../ordering/ordering.hpp"
#include "../seqbuilder/seqbuilder.hpp"
#include "../sais/sais.hpp"
#include "../fmindex/fmindex.hpp"
#include "../streamencode/streamencode.hpp"
#include "../windows/windows.hpp"
#include "../bitvectors/bitvectors.hpp"
#include "../sarange/sarange.hpp"
#include "../seal/seal.hpp"

#include "bamsix/config.hpp"

#include <chrono>
#include <cstring>
#include <iostream>

namespace bamsix {

void BuildIndex(const std::vector<std::string>& bam_paths,
                const std::string& output_path,
                const BuildConfig& config) {
    auto t0 = std::chrono::steady_clock::now();
    auto elapsed = [&]() {
        auto now = std::chrono::steady_clock::now();
        return std::chrono::duration_cast<std::chrono::milliseconds>(now - t0).count();
    };

    // ─── Stage 1: Ingestion ──────────────────────────────────────────────────
    std::cerr << "[build] Stage 1: Ingesting BAM files...\n";
    auto ingest = IngestBams(bam_paths);
    std::cerr << "[build]   " << ingest.reads.size() << " reads ingested from "
              << bam_paths.size() << " file(s) (" << elapsed() << " ms)\n";

    if (ingest.reads.empty()) {
        throw Error{ErrorCode::INVALID_BAM_INPUT,
                    "No reads passed the inclusion rule"};
    }

    // H1 fix: reject --enable-bidirectional in v1.0 (reverse FM-index not implemented)
    if (config.enable_bidirectional) {
        throw Error{ErrorCode::BUILD_VALIDATION_FAILED,
                    "--enable-bidirectional is not implemented in v1.0. "
                    "Reverse FM-index build/serialize/deserialize requires v2.0"};
    }

    // ─── Stage 2: Ordering ───────────────────────────────────────────────────
    std::cerr << "[build] Stage 2: Ordering reads...\n";
    auto ordering = OrderReads(ingest);
    std::cerr << "[build]   " << ordering.reads.size() << " reads ordered (" << elapsed() << " ms)\n";

    // ─── Stage 3: Sequence construction ──────────────────────────────────────
    std::cerr << "[build] Stage 3: Building concatenated sequence S...\n";
    auto bundle = BuildSequence(ordering.reads);
    std::cerr << "[build]   |S| = " << bundle.S.size()
              << ", N = " << bundle.readStarts.size()
              << " (" << elapsed() << " ms)\n";

    // ─── Stage 4: SA-IS construction ─────────────────────────────────────────
    std::cerr << "[build] Stage 4: Computing suffix array (SA-IS)...\n";
    auto sais = ComputeSuffixArray(bundle);
    std::cerr << "[build]   |SA| = " << sais.SA.size()
              << ", sentinel_row = " << sais.sentinel_row
              << " (" << elapsed() << " ms)\n";

    // ─── Stage 4b: ISA samples (optional, Architecture §4.4 step 5) ─────────
    if (config.sample_step_s_prime > 0) {
        std::cerr << "[build] Stage 4b: Computing ISA samples (s'="
                  << config.sample_step_s_prime << ")...\n";
        ComputeISASamples(sais, config.sample_step_s_prime);
        std::cerr << "[build]   ISA samples: " << sais.ISA_samples.size()
                  << " (" << elapsed() << " ms)\n";
    }

    // ─── Stage 5b: FM-index construction (before S_seq — need stripped BWT) ──
    std::cerr << "[build] Stage 5b: Building FM-index...\n";
    FMIndexEngine fm;
    fm.Build(sais.BWT, sais.SA, sais.sentinel_row, config.sample_step_s,
             bundle.S.size());
    if (!sais.ISA_samples.empty()) {
        fm.SetISASamples(std::move(sais.ISA_samples), config.sample_step_s_prime);
    }
    std::cerr << "[build]   SA samples: " << fm.SASamples().size()
              << (fm.HasISASamples() ? ", ISA samples loaded" : "")
              << " (" << elapsed() << " ms)\n";

    // OOM FIX: sais.SA is ~16.8 GB for a 2B text. Now that we have sa_samples,
    // we must destroy the full SA immediately to free RAM before Stage 5a.
    std::vector<int64_t>().swap(sais.SA);
    std::vector<uint8_t>().swap(sais.BWT);

    // ─── Stage 5a: S_seq encoding ────────────────────────────────────────────
    // C1 fix: encode the sentinel-stripped BWT (|S| entries, codes 0–5 only).
    // The full BWT contains CODE_SENT=6 which is outside MTF alphabet SIGMA=6.
    // fm.SerializeBWT() returns BWT with sentinel row removed.
    std::cerr << "[build] Stage 5a: Encoding S_seq...\n";
    auto stripped_bwt = fm.SerializeBWT();
    auto seq_encoded = EncodeSeqStream(stripped_bwt, config.entropy_order_k);
    std::vector<uint8_t>().swap(stripped_bwt);
    std::cerr << "[build]   S_seq payload: " << seq_encoded.payload.size() << " bytes"
              << " (" << elapsed() << " ms)\n";

    // ─── Stage 6: Stream encoding (S_qual, S_meta, S_map) ───────────────────
    std::cerr << "[build] Stage 6: Encoding auxiliary streams...\n";
    auto qual_encoded = EncodeQualStream(ordering.reads, config.qual_codec,
                                          config.qual_lossy_bins,
                                          config.qual_block_size);
    auto meta_encoded = EncodeMetaStream(ordering.reads, config.meta_codec);
    auto map_encoded  = EncodeMapStream(ordering.reads, config.map_codec);
    std::cerr << "[build]   S_qual: " << qual_encoded.payload.size() << " bytes\n";
    std::cerr << "[build]   S_meta: " << meta_encoded.payload.size() << " bytes\n";
    std::cerr << "[build]   S_map:  " << map_encoded.payload.size() << " bytes"
              << " (" << elapsed() << " ms)\n";

    // ─── Stage 7: Window construction ────────────────────────────────────────
    std::cerr << "[build] Stage 7: Building window table...\n";
    auto windows = BuildWindows(ordering.reads, bundle, config.window_size_T);
    std::cerr << "[build]   " << windows.size() << " windows (" << elapsed() << " ms)\n";

    // ─── Stage 8: Bitvector construction ─────────────────────────────────────
    std::cerr << "[build] Stage 8: Building bitvectors...\n";
    auto bv = BuildBitvectors(bundle, windows);
    std::cerr << "[build]   B_read popcount = " << bv.B_read.PopCount()
              << ", B_window popcount = " << bv.B_window.PopCount()
              << " (" << elapsed() << " ms)\n";

    // ─── Stage 5c: SARange construction (optional, ENHANCED tier) ────────────
    SARange sarange;
    if (config.enable_sarange) {
        std::cerr << "[build] Stage 5c: Building SARange wavelet tree (ENHANCED)...\n";
        sarange.Build(fm.SASamples(), bundle.S.size(), config.sample_step_s);
        std::cerr << "[build]   SARange built over " << fm.SASamples().size()
                  << " SA samples (" << elapsed() << " ms)\n";
    }

    // ─── Stage 9-10: Validation + Sealing ────────────────────────────────────
    std::cerr << "[build] Stage 9-10: Sealing .bsi file...\n";

    // Populate header
    BsiHeader header;
    std::memcpy(header.magic, "BMSI", 4);
    header.version = 6;
    std::strncpy(header.bamsix_version, BAMSIX_VERSION, 15);

    // Provenance fields per Contract §9.2
    header.build_timestamp_utc = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::seconds>(
            std::chrono::system_clock::now().time_since_epoch()).count());
#if defined(__linux__)
    header.host_os_id = 1;  // Linux
#elif defined(__APPLE__)
    header.host_os_id = 2;  // macOS
#elif defined(_WIN32)
    header.host_os_id = 3;  // Windows
#else
    header.host_os_id = 0;  // Unknown
#endif
#if defined(__x86_64__) || defined(_M_X64)
    header.cpu_arch_id = 1;  // x86_64
#elif defined(__aarch64__) || defined(_M_ARM64)
    header.cpu_arch_id = 2;  // ARM64
#else
    header.cpu_arch_id = 0;  // Unknown
#endif

    header.is_lossless = (config.qual_lossy_bins == 0) ? 1 : 0;
    header.source_file_count = ingest.source_file_count;
    header.source_manifest_hash = ingest.source_manifest_hash;
    header.ordering_hash = ordering.ordering_hash;
    header.S_length = bundle.S.size();
    header.N_reads = ordering.reads.size();
    header.N_windows = static_cast<uint32_t>(windows.size());
    header.sample_step_s = static_cast<uint32_t>(config.sample_step_s);
    header.has_isa_samples = (config.sample_step_s_prime > 0) ? 1 : 0;
    header.sample_step_s_prime = static_cast<uint32_t>(config.sample_step_s_prime);
    header.enable_sarange = config.enable_sarange ? 1 : 0;
    header.shared_bwt = config.shared_bwt ? 1 : 0;
    header.enable_bidirectional = config.enable_bidirectional ? 1 : 0;
    header.recommended_seed_length = config.recommended_seed_length;
    header.window_size_T = config.window_size_T;
    header.entropy_order_k = config.entropy_order_k;
    header.qual_codec_id = static_cast<uint8_t>(config.qual_codec);
    header.qual_lossy_bins = config.qual_lossy_bins;
    header.meta_codec_id = static_cast<uint8_t>(config.meta_codec);
    header.map_codec_id = static_cast<uint8_t>(config.map_codec);
    header.strand_mode = static_cast<uint8_t>(config.strand_mode);
    header.sentinel_row = fm.SentinelRow();
    header.chrom_count = static_cast<uint32_t>(ingest.chrom_names.size());
    for (uint32_t i = 0; i < ingest.chrom_names.size(); ++i) {
        header.chrom_name_table.push_back({i, ingest.chrom_names[i]});
    }
    header.seq_block_size = config.seq_block_size;
    header.qual_block_size = config.qual_block_size;
    header.allow_parallel_sa = config.allow_parallel_sa ? 1 : 0;
    header.reference_based_encoding = config.reference_fasta_path.empty() ? 0 : 1;
    header.reference_sha256 = config.reference_sha256;

    SealInput seal_input;
    seal_input.header = header;
    seal_input.seq = std::move(seq_encoded);
    seal_input.qual = std::move(qual_encoded);
    seal_input.meta = std::move(meta_encoded);
    seal_input.map = std::move(map_encoded);
    seal_input.bv = std::move(bv);
    seal_input.windows = std::move(windows);
    seal_input.reads = std::move(ordering.reads);

    // V5 ENHANCED features
    if (fm.HasISASamples()) {
        seal_input.isa_samples = fm.ISASamples(); // 264 MB copy is fine
        seal_input.isa_step = fm.ISAStep();
    }
    if (sarange.IsBuilt()) {
        seal_input.sarange = std::move(sarange);
    }
    seal_input.fm = std::move(fm); // Move FM last since we extracted ISA samples from it

    WriteBsi(output_path, seal_input);

    auto total_ms = elapsed();
    std::cerr << "[build] Done. Output: " << output_path
              << " (" << total_ms << " ms total)\n";
}

}  // namespace bamsix
