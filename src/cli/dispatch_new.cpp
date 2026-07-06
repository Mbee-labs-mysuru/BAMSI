#include "bamsix/cli/dispatch.hpp"
#include "build.hpp"
#include "../query/query.hpp"
#include "../fmindex/fmindex.hpp"
#include "../bitvectors/bitvectors.hpp"
#include "../format/format.hpp"
#include "../streamencode/streamencode.hpp"
#include "bamsix/config.hpp"

#include <htslib/sam.h>
#include <htslib/hts.h>

#include <openssl/sha.h>

#include <iostream>
#include <cstring>
#include <string>
#include <unordered_set>
#include <vector>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <fstream>
#include <chrono>

namespace bamsix {

namespace {

// Forward declarations for utility parsers (defined later, used by early handlers)
std::string ParseStringArg(int argc, char** argv, const std::string& flag);
uint64_t ParseU64Arg(int argc, char** argv, const std::string& flag, uint64_t def);

/// Parse --input values (can appear multiple times) and positional BAM args.
std::vector<std::string> ParseInputs(int argc, char** argv) {
    // All flags that take a value argument (must be skipped when parsing positionals)
    static const std::unordered_set<std::string> kValueFlags = {
        "--output", "-o", "--window-size", "--sample-step", "--entropy-k",
        "--lossy-bins", "--isa-step", "--seq-block-size", "--qual-block-size",
        "--seed-length", "--strand", "--reference", "--index", "--pattern",
        "--chrom", "--start", "--end", "--threshold", "--region", "--streams",
        "--read-ids", "--read-id", "--qual-codec", "--meta-codec", "--map-codec",
        "--threads", "--bed", "--overlap",
    };
    std::vector<std::string> inputs;
    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--input" && i + 1 < argc) {
            inputs.push_back(argv[++i]);
        } else if (kValueFlags.count(arg)) {
            ++i;  // skip the value
        } else if (arg[0] == '-') {
            // skip boolean flags
        } else {
            inputs.push_back(arg);
        }
    }
    return inputs;
}

/// Parse --output or -o value
std::string ParseOutput(int argc, char** argv, const std::string& def = "") {
    for (int i = 0; i < argc; ++i) {
        std::string a = argv[i];
        if ((a == "--output" || a == "-o") && i + 1 < argc) {
            return argv[i + 1];
        }
    }
    return def;
}

/// Parse --index value (also accepts positional .bsi arg)
std::string ParseIndex(int argc, char** argv) {
    for (int i = 2; i < argc; ++i) {
        if (std::string(argv[i]) == "--index" && i + 1 < argc) {
            return argv[i + 1];
        }
    }
    // Fallback: first positional arg ending in .bsi
    for (int i = 2; i < argc; ++i) {
        std::string a = argv[i];
        if (a[0] != '-' && a.size() > 4 &&
            a.substr(a.size() - 4) == ".bsi") {
            return a;
        }
    }
    return "";
}

/// Parse --pattern value (string to code vector)
std::vector<uint8_t> ParsePattern(int argc, char** argv) {
    std::string pat;
    for (int i = 0; i < argc; ++i) {
        if (std::string(argv[i]) == "--pattern" && i + 1 < argc) {
            pat = argv[i + 1];
            break;
        }
    }
    std::vector<uint8_t> codes;
    for (char c : pat) {
        switch (c) {
            case 'A': case 'a': codes.push_back(CODE_A); break;
            case 'C': case 'c': codes.push_back(CODE_C); break;
            case 'G': case 'g': codes.push_back(CODE_G); break;
            case 'T': case 't': codes.push_back(CODE_T); break;
            default:            codes.push_back(CODE_N); break;
        }
    }
    return codes;
}

bool HasFlag(int argc, char** argv, const std::string& flag) {
    for (int i = 0; i < argc; ++i) {
        if (std::string(argv[i]) == flag) return true;
    }
    return false;
}

// ─── Build ──────────────────────────────────────────────────────────────────

Status HandleBuild(int argc, char** argv) {
    auto inputs = ParseInputs(argc, argv);
    auto output = ParseOutput(argc, argv, "output.bsi");

    if (inputs.empty()) {
        std::cerr << "error: --input required\n";
        return Status(StatusCode::kInvalidArgument, "Missing --input");
    }

    BuildConfig config;
    for (int i = 0; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--window-size" && i + 1 < argc)
            config.window_size_T = std::stoull(argv[++i]);
        else if (arg == "--sample-step" && i + 1 < argc)
            config.sample_step_s = std::stoull(argv[++i]);
        else if (arg == "--entropy-k" && i + 1 < argc)
            config.entropy_order_k = static_cast<uint8_t>(std::stoul(argv[++i]));
        else if (arg == "--lossy-bins" && i + 1 < argc)
            config.qual_lossy_bins = static_cast<uint8_t>(std::stoul(argv[++i]));
        else if (arg == "--isa-step" && i + 1 < argc)
            config.sample_step_s_prime = std::stoull(argv[++i]);
        else if (arg == "--seq-block-size" && i + 1 < argc)
            config.seq_block_size = static_cast<uint32_t>(std::stoul(argv[++i]));
        else if (arg == "--qual-block-size" && i + 1 < argc)
            config.qual_block_size = static_cast<uint32_t>(std::stoul(argv[++i]));
        else if (arg == "--seed-length" && i + 1 < argc)
            config.recommended_seed_length = static_cast<uint8_t>(std::stoul(argv[++i]));
        else if (arg == "--strand" && i + 1 < argc) {
            std::string s = argv[++i];
            if (s == "single") config.strand_mode = StrandMode::SingleStrand;
            else config.strand_mode = StrandMode::StrandComplete;
        }
        else if (arg == "--reference" && i + 1 < argc)
            config.reference_fasta_path = argv[++i];
        else if (arg == "--lossless")
            config.qual_lossy_bins = 0;
        else if (arg == "--lossy")
            config.qual_lossy_bins = config.qual_lossy_bins > 0 ? config.qual_lossy_bins : 8;
        else if (arg == "--enable-sarange")
            config.enable_sarange = true;
        else if (arg == "--enable-bidirectional")
            config.enable_bidirectional = true;
        else if (arg == "--shared-bwt")
            config.shared_bwt = true;
        else if (arg == "--parallel-sa")
            config.allow_parallel_sa = true;
        else if (arg == "--qual-codec" && i + 1 < argc) {
            std::string codec = argv[++i];
            if (codec == "RANGE_CYCLE")       config.qual_codec = QualCodec::RANGE_CYCLE;
            else if (codec == "RANS_DELTA")   config.qual_codec = QualCodec::RANS_DELTA;
            else if (codec == "ZSTD_DICT")    config.qual_codec = QualCodec::ZSTD_DICT;
            else if (codec == "BINNED_RANGE") config.qual_codec = QualCodec::BINNED_RANGE;
            else { std::cerr << "error: unknown --qual-codec: " << codec << "\n"; }
        }
        else if (arg == "--meta-codec" && i + 1 < argc) {
            std::string codec = argv[++i];
            if (codec == "TYPED_SPLIT")       config.meta_codec = MetaCodec::TYPED_SPLIT;
            else if (codec == "ZSTD_FALLBACK") config.meta_codec = MetaCodec::ZSTD_FALLBACK;
            else { std::cerr << "error: unknown --meta-codec: " << codec << "\n"; }
        }
        else if (arg == "--map-codec" && i + 1 < argc) {
            std::string codec = argv[++i];
            if (codec == "DELTA_RANGE")       config.map_codec = MapCodec::DELTA_RANGE;
            else if (codec == "RAW")          config.map_codec = MapCodec::RAW;
            else { std::cerr << "error: unknown --map-codec: " << codec << "\n"; }
        }
        else if (arg == "--threads" && i + 1 < argc) {
            config.build_threads = static_cast<uint32_t>(std::stoul(argv[++i]));
        }
    }

    try {
        BuildIndex(inputs, output, config);
        return Status::Ok();
    } catch (const Error& e) {
        std::cerr << "error: " << e.message << "\n";
        return Status(StatusCode::kInvalidArgument, e.message);
    }
}

// ─── Count ──────────────────────────────────────────────────────────────────

Status HandleCount(int argc, char** argv) {
    auto index_path = ParseIndex(argc, argv);
    auto pattern = ParsePattern(argc, argv);
    bool benchmark = HasFlag(argc, argv, "--benchmark");

    if (index_path.empty()) {
        std::cerr << "error: --index <file.bsi> required\n";
        return Status(StatusCode::kInvalidArgument, "Missing --index");
    }
    if (pattern.empty()) {
        std::cerr << "error: --pattern <ACGT...> required\n";
        return Status(StatusCode::kInvalidArgument, "Missing --pattern");
    }

    try {
        // Phase 1: Index loading (I/O-dominated)
        auto t0 = std::chrono::steady_clock::now();
        auto idx = ReadBsi(index_path);
        auto t1 = std::chrono::steady_clock::now();

        StrandMode mode = static_cast<StrandMode>(idx.header.strand_mode);
        // --strand override
        auto strand_arg = ParseStringArg(argc, argv, "--strand");
        if (strand_arg == "single" || strand_arg == "SingleStrand")
            mode = StrandMode::SingleStrand;
        else if (strand_arg == "complete" || strand_arg == "StrandComplete")
            mode = StrandMode::StrandComplete;

        // --overlap override (default: overlapping; "none" = non-overlapping)
        OverlapMode overlap = OverlapMode::All;
        auto overlap_arg = ParseStringArg(argc, argv, "--overlap");
        if (overlap_arg == "none" || overlap_arg == "nonoverlapping")
            overlap = OverlapMode::None;

        // Phase 2: Backward search (O(|P|) algorithmic phase)
        auto t2 = std::chrono::steady_clock::now();
        uint64_t count = GlobalCount(pattern, idx.fm, mode, overlap);
        auto t3 = std::chrono::steady_clock::now();

        double load_s = std::chrono::duration<double>(t1 - t0).count();
        double search_s = std::chrono::duration<double>(t3 - t2).count();
        double total_s = std::chrono::duration<double>(t3 - t0).count();

        if (HasFlag(argc, argv, "--json")) {
            std::cout << "{\"count\":" << count;
            if (benchmark) {
                std::cout << ",\"load_time_s\":" << std::fixed << std::setprecision(6) << load_s
                          << ",\"search_time_s\":" << search_s
                          << ",\"total_time_s\":" << total_s;
            }
            std::cout << "}\n";
        } else {
            std::cout << count << "\n";
            if (benchmark) {
                std::cerr << "[benchmark] load=" << std::fixed << std::setprecision(6) << load_s
                          << "s search=" << search_s
                          << "s total=" << total_s << "s\n";
            }
        }
        return Status::Ok();
    } catch (const Error& e) {
        std::cerr << "error: " << e.message << "\n";
        return Status(StatusCode::kInvalidArgument, e.message);
    }
}

// ─── Exists ─────────────────────────────────────────────────────────────────

Status HandleExists(int argc, char** argv) {
    auto index_path = ParseIndex(argc, argv);
    auto pattern = ParsePattern(argc, argv);

    if (index_path.empty() || pattern.empty()) {
        std::cerr << "error: --index and --pattern required\n";
        return Status(StatusCode::kInvalidArgument, "Missing args");
    }

    try {
        auto idx = ReadBsi(index_path);
        StrandMode mode = static_cast<StrandMode>(idx.header.strand_mode);
        // --strand override
        auto strand_arg = ParseStringArg(argc, argv, "--strand");
        if (strand_arg == "single" || strand_arg == "SingleStrand")
            mode = StrandMode::SingleStrand;
        else if (strand_arg == "complete" || strand_arg == "StrandComplete")
            mode = StrandMode::StrandComplete;
        bool exists = GlobalExists(pattern, idx.fm, mode);

        if (HasFlag(argc, argv, "--json")) {
            std::cout << "{\"exists\":" << (exists ? "true" : "false") << "}\n";
        } else {
            std::cout << (exists ? "true" : "false") << "\n";
        }
        return Status::Ok();
    } catch (const Error& e) {
        std::cerr << "error: " << e.message << "\n";
        return Status(StatusCode::kInvalidArgument, e.message);
    }
}

// ─── Locate ─────────────────────────────────────────────────────────────────

Status HandleLocate(int argc, char** argv) {
    auto index_path = ParseIndex(argc, argv);
    auto pattern = ParsePattern(argc, argv);
    // Contract §10.2: --sort-output (also accept --sorted for compat)
    bool sorted = HasFlag(argc, argv, "--sort-output") || HasFlag(argc, argv, "--sorted");
    bool bed_fmt = HasFlag(argc, argv, "--bed");
    auto output_path = ParseOutput(argc, argv, "");

    if (index_path.empty() || pattern.empty()) {
        std::cerr << "error: --index and --pattern required\n";
        return Status(StatusCode::kInvalidArgument, "Missing args");
    }

    try {
        auto idx = ReadBsi(index_path);
        StrandMode mode = static_cast<StrandMode>(idx.header.strand_mode);
        // --strand override
        auto strand_arg = ParseStringArg(argc, argv, "--strand");
        if (strand_arg == "single" || strand_arg == "SingleStrand")
            mode = StrandMode::SingleStrand;
        else if (strand_arg == "complete" || strand_arg == "StrandComplete")
            mode = StrandMode::StrandComplete;

        // Use lazy decode from streams (Contract §2.3 Table B — production path)
        // Falls back to legacy in-memory path if directories are empty
        std::vector<Match> matches;
        if (!idx.meta_directory.empty() && !idx.map_directory.empty()) {
            // Build dense absolute-position table ONCE for O(1) delta resolution
            // instead of O(N) per-occurrence backward+forward scan
            auto full_abs_table = BuildFullAbsPosTable(
                idx.map_payload, idx.map_directory, idx.header.map_codec_id);
            matches = LocateLazy(pattern, idx.fm, idx.bv.B_read,
                                 idx.meta_payload, idx.map_payload,
                                 idx.meta_directory, idx.map_directory,
                                 idx.header.meta_codec_id, idx.header.map_codec_id,
                                 idx.chrom_names, mode, sorted, &full_abs_table);
        } else {
            matches = Locate(pattern, idx.fm, idx.bv.B_read,
                             idx.reads, idx.chrom_names, mode, sorted);
        }

        // Output destination
        std::ofstream ofs;
        std::ostream* out = &std::cout;
        if (!output_path.empty()) {
            ofs.open(output_path);
            if (!ofs.is_open()) {
                std::cerr << "error: cannot open output file: " << output_path << "\n";
                return Status(StatusCode::kInvalidArgument, "Cannot open output");
            }
            out = &ofs;
        }

        if (HasFlag(argc, argv, "--json")) {
            *out << "{\"mode\":\"" << (sorted ? "sorted" : "streaming")
                 << "\",\"matches\":[";
            for (size_t i = 0; i < matches.size(); ++i) {
                if (i > 0) *out << ",";
                const auto& m = matches[i];
                *out << "{\"strand\":\"" << (m.query_strand == QueryStrand::Forward ? "+" : "-")
                     << "\",\"chrom\":\"" << m.chrom
                     << "\",\"p_min\":" << m.p_min
                     << ",\"p_max\":" << m.p_max
                     << ",\"read_id\":" << m.read_id << "}";
            }
            *out << "],\"count\":" << matches.size() << "}\n";
        } else if (bed_fmt) {
            // BED format: chrom  start(0-based)  end(exclusive)  name  score  strand
            for (const auto& m : matches) {
                *out << m.chrom << "\t" << (m.p_min - 1) << "\t" << m.p_max
                     << "\tread_" << m.read_id << "\t0\t"
                     << (m.query_strand == QueryStrand::Forward ? "+" : "-") << "\n";
            }
        } else {
            // TSV header
            *out << "strand\tchrom\tp_min\tp_max\tread_id\n";
            for (const auto& m : matches) {
                *out << (m.query_strand == QueryStrand::Forward ? "+" : "-")
                     << "\t" << m.chrom
                     << "\t" << m.p_min
                     << "\t" << m.p_max
                     << "\t" << m.read_id << "\n";
            }
        }
        return Status::Ok();
    } catch (const Error& e) {
        std::cerr << "error: " << e.message << "\n";
        return Status(StatusCode::kInvalidArgument, e.message);
    }
}

// ─── Verify ─────────────────────────────────────────────────────────────────

Status HandleVerify(int argc, char** argv) {
    auto index_path = ParseIndex(argc, argv);
    bool strict = HasFlag(argc, argv, "--strict");

    if (index_path.empty()) {
        std::cerr << "error: --index <file.bsi> required\n";
        return Status(StatusCode::kInvalidArgument, "Missing --index");
    }

    try {
        bool ok = VerifyBsi(index_path);
        std::string detail;

        if (strict && ok) {
            // --strict: also open index and verify all section checksums
            // plus ordering_hash and source_manifest_hash are present
            auto idx = ReadBsi(index_path);
            const auto& h = idx.header;

            // Check ordering_hash is non-zero
            bool hash_nonzero = false;
            for (auto b : h.ordering_hash) {
                if (b != 0) { hash_nonzero = true; break; }
            }
            if (!hash_nonzero) {
                detail += "  WARNING: ordering_hash is all-zero\n";
            }

            // Check source_manifest_hash is non-zero
            hash_nonzero = false;
            for (auto b : h.source_manifest_hash) {
                if (b != 0) { hash_nonzero = true; break; }
            }
            if (!hash_nonzero) {
                detail += "  WARNING: source_manifest_hash is all-zero\n";
            }

            // Verify read count consistency
            if (idx.reads.size() != h.N_reads) {
                detail += "  FAIL: N_reads mismatch (header=" + std::to_string(h.N_reads)
                        + " actual=" + std::to_string(idx.reads.size()) + ")\n";
                ok = false;
            }

            // Verify window count
            if (idx.windows.size() != h.N_windows) {
                detail += "  FAIL: N_windows mismatch (header=" + std::to_string(h.N_windows)
                        + " actual=" + std::to_string(idx.windows.size()) + ")\n";
                ok = false;
            }

            // Verify chrom table
            if (idx.chrom_names.size() != h.chrom_count) {
                detail += "  FAIL: chrom_count mismatch\n";
                ok = false;
            }

            // Re-compute ordering_hash from stored reads (Architecture §4.2)
            // ordering_hash = SHA-256(for i=0..N-1: uint32_le(chrom_id) || uint64_le(pos)
            //                         || uint32_le(source_file_id) || uint64_le(bam_offset))
            {
                SHA256_CTX sha_ctx;
                SHA256_Init(&sha_ctx);
                for (size_t i = 0; i < idx.reads.size(); ++i) {
                    const auto& rd = idx.reads[i];
                    uint32_t cid = rd.chrom_id;
                    uint64_t pos = rd.pos;
                    uint32_t sfid = rd.source_file_id;
                    uint64_t boff = rd.bam_offset;
                    SHA256_Update(&sha_ctx, &cid, 4);
                    SHA256_Update(&sha_ctx, &pos, 8);
                    SHA256_Update(&sha_ctx, &sfid, 4);
                    SHA256_Update(&sha_ctx, &boff, 8);
                }
                std::array<uint8_t, 32> recomputed;
                SHA256_Final(recomputed.data(), &sha_ctx);

                if (recomputed != h.ordering_hash) {
                    detail += "  FAIL: ordering_hash mismatch (stored reads inconsistent with header)\n";
                    ok = false;
                } else {
                    detail += "  OK: ordering_hash re-verified from stored reads\n";
                }
            }
        }

        if (HasFlag(argc, argv, "--json")) {
            std::cout << "{\"valid\":" << (ok ? "true" : "false")
                      << ",\"strict\":" << (strict ? "true" : "false") << "}\n";
        } else {
            if (ok) {
                std::cout << "PASS: " << index_path << " integrity verified"
                          << (strict ? " (strict)" : "") << "\n";
            } else {
                std::cout << "FAIL: " << index_path << " verification failed\n";
            }
            if (!detail.empty()) {
                std::cout << detail;
            }
        }
        return ok ? Status::Ok()
                  : Status(StatusCode::kInvalidArgument, "Verification failed");
    } catch (const Error& e) {
        std::cerr << "error: " << e.message << "\n";
        return Status(StatusCode::kInvalidArgument, e.message);
    }
}

// ─── Info ───────────────────────────────────────────────────────────────────

Status HandleInfo(int argc, char** argv) {
    auto index_path = ParseIndex(argc, argv);

    if (index_path.empty()) {
        std::cerr << "error: --index <file.bsi> required\n";
        return Status(StatusCode::kInvalidArgument, "Missing --index");
    }

    try {
        auto idx = ReadBsi(index_path);
        const auto& h = idx.header;
        bool json = HasFlag(argc, argv, "--json");

        // Compute BWT run count: r = number of maximal equal-letter runs
        // A single linear scan over the stored BWT — O(|S|) time.
        uint64_t bwt_total_len = idx.fm.SLen() + 1;  // BWT has |S|+1 entries
        uint64_t bwt_runs = 1;  // first character starts a run
        for (uint64_t i = 1; i < bwt_total_len; ++i) {
            if (idx.fm.BwtAt(i) != idx.fm.BwtAt(i - 1))
                ++bwt_runs;
        }
        double s_over_r = (bwt_runs > 0) ? (double)idx.fm.SLen() / bwt_runs : 0.0;

        // Ensure bamsix_version is null-terminated (Mo7 safety)
        char safe_version[17] = {};
        std::memcpy(safe_version, h.bamsix_version, 16);
        safe_version[16] = '\0';

        // Helper: format 32-byte hash as hex
        auto hex32 = [](const std::array<uint8_t, 32>& hash) -> std::string {
            std::ostringstream oss;
            for (auto b : hash) oss << std::hex << std::setfill('0') << std::setw(2) << (int)b;
            return oss.str();
        };

        if (json) {
            // Full JSON output per Contract §9.2 — all 30+ fields
            std::cout << "{\n";
            std::cout << "  \"magic\": \"BMSI\",\n";
            std::cout << "  \"format_version\": " << h.version << ",\n";
            std::cout << "  \"bamsix_version\": \"" << safe_version << "\",\n";
            std::cout << "  \"host_os_id\": " << (int)h.host_os_id << ",\n";
            std::cout << "  \"cpu_arch_id\": " << (int)h.cpu_arch_id << ",\n";
            std::cout << "  \"build_timestamp_utc\": " << h.build_timestamp_utc << ",\n";
            std::cout << "  \"is_lossless\": " << (h.is_lossless ? "true" : "false") << ",\n";
            std::cout << "  \"source_file_count\": " << h.source_file_count << ",\n";
            std::cout << "  \"source_manifest_hash\": \"" << hex32(h.source_manifest_hash) << "\",\n";
            std::cout << "  \"ordering_hash\": \"" << hex32(h.ordering_hash) << "\",\n";
            std::cout << "  \"S_length\": " << h.S_length << ",\n";
            std::cout << "  \"N_reads\": " << h.N_reads << ",\n";
            std::cout << "  \"N_windows\": " << h.N_windows << ",\n";
            std::cout << "  \"sample_step_s\": " << h.sample_step_s << ",\n";
            std::cout << "  \"has_isa_samples\": " << (h.has_isa_samples ? "true" : "false") << ",\n";
            std::cout << "  \"sample_step_s_prime\": " << h.sample_step_s_prime << ",\n";
            std::cout << "  \"enable_sarange\": " << (h.enable_sarange ? "true" : "false") << ",\n";
            std::cout << "  \"sarange_variant\": " << (int)h.sarange_variant << ",\n";
            std::cout << "  \"shared_bwt\": " << (h.shared_bwt ? "true" : "false") << ",\n";
            std::cout << "  \"enable_bidirectional\": " << (h.enable_bidirectional ? "true" : "false") << ",\n";
            std::cout << "  \"recommended_seed_length\": " << (int)h.recommended_seed_length << ",\n";
            std::cout << "  \"window_size_T\": " << h.window_size_T << ",\n";
            std::cout << "  \"entropy_order_k\": " << (int)h.entropy_order_k << ",\n";
            std::cout << "  \"qual_codec_id\": " << (int)h.qual_codec_id << ",\n";
            std::cout << "  \"qual_lossy_bins\": " << (int)h.qual_lossy_bins << ",\n";
            std::cout << "  \"meta_codec_id\": " << (int)h.meta_codec_id << ",\n";
            std::cout << "  \"map_codec_id\": " << (int)h.map_codec_id << ",\n";
            std::cout << "  \"strand_mode\": " << (int)h.strand_mode << ",\n";
            std::cout << "  \"sentinel_row\": " << h.sentinel_row << ",\n";
            std::cout << "  \"seq_block_size\": " << h.seq_block_size << ",\n";
            std::cout << "  \"qual_block_size\": " << h.qual_block_size << ",\n";
            std::cout << "  \"allow_parallel_sa\": " << (h.allow_parallel_sa ? "true" : "false") << ",\n";
            std::cout << "  \"reference_based_encoding\": " << (h.reference_based_encoding ? "true" : "false") << ",\n";
            if (h.reference_based_encoding) {
                std::cout << "  \"reference_sha256\": \"" << hex32(h.reference_sha256) << "\",\n";
            }
            std::cout << "  \"sa_samples_count\": " << idx.fm.SASamples().size() << ",\n";
            std::cout << "  \"chrom_count\": " << h.chrom_count << ",\n";
            std::cout << "  \"chrom_name_table\": [";
            for (uint32_t i = 0; i < h.chrom_count; ++i) {
                if (i > 0) std::cout << ", ";
                std::cout << "{\"chrom_id\": " << h.chrom_name_table[i].chrom_id
                          << ", \"name\": \"" << h.chrom_name_table[i].name << "\"}";
            }
            std::cout << "],\n";
            // M6 fix: add mode (BASE/ENHANCED) and human-readable strand_mode_name
            std::string mode_str = (h.enable_sarange || h.enable_bidirectional) ? "ENHANCED" : "BASE";
            std::string strand_name = (h.strand_mode == 0) ? "strand_complete" : "single_strand";
            std::cout << "  \"mode\": \"" << mode_str << "\",\n";
            std::cout << "  \"strand_mode_name\": \"" << strand_name << "\",\n";
            std::cout << "  \"bwt_runs\": " << bwt_runs << ",\n";
            std::cout << "  \"S_over_r\": " << std::fixed << std::setprecision(2) << s_over_r << "\n";
            std::cout << "}\n";
        } else {
            std::cout << "BAMSIX Index: " << index_path << "\n";
            std::cout << "  Format version:      " << h.version << "\n";
            std::cout << "  BAMSIX version:       " << safe_version << "\n";
            std::cout << "  Lossless:            " << (h.is_lossless ? "yes" : "no") << "\n";
            std::cout << "  |S|:                 " << h.S_length << "\n";
            std::cout << "  Reads:               " << h.N_reads << "\n";
            std::cout << "  Windows:             " << h.N_windows << "\n";
            std::cout << "  SA sample step:      " << h.sample_step_s << "\n";
            std::cout << "  Window size T:       " << h.window_size_T << "\n";
            std::cout << "  Entropy order k:     " << (int)h.entropy_order_k << "\n";
            std::cout << "  Shared BWT:          " << (h.shared_bwt ? "yes" : "no") << "\n";
            std::cout << "  Strand mode:         " << (h.strand_mode == 0 ? "strand-complete" : "single-strand") << "\n";
            std::cout << "  SARange (ENHANCED):  " << (h.enable_sarange ? "yes" : "no") << "\n";
            std::cout << "  Bidirectional FM:    " << (h.enable_bidirectional ? "yes" : "no") << "\n";
            std::cout << "  ISA samples:         " << (h.has_isa_samples ? "yes" : "no") << "\n";
            std::cout << "  Qual codec:          0x" << std::hex << (int)h.qual_codec_id << std::dec << "\n";
            std::cout << "  Meta codec:          0x" << std::hex << (int)h.meta_codec_id << std::dec << "\n";
            std::cout << "  Map codec:           0x" << std::hex << (int)h.map_codec_id << std::dec << "\n";
            if (h.qual_lossy_bins > 0) {
                std::cout << "  Lossy bins:          " << (int)h.qual_lossy_bins << "\n";
            }
            std::cout << "  Source manifest:     " << hex32(h.source_manifest_hash) << "\n";
            std::cout << "  Ordering hash:       " << hex32(h.ordering_hash) << "\n";
            std::cout << "  Chromosomes:         " << h.chrom_count << "\n";
            for (const auto& ce : h.chrom_name_table) {
                std::cout << "    [" << ce.chrom_id << "] " << ce.name << "\n";
            }
            std::cout << "  SA samples:          " << idx.fm.SASamples().size() << "\n";
            std::cout << "  Sentinel row:        " << h.sentinel_row << "\n";
            std::cout << "  BWT runs (r):        " << bwt_runs << "\n";
            std::cout << "  |S|/r ratio:         " << std::fixed << std::setprecision(2) << s_over_r << "\n";
        }
        return Status::Ok();
    } catch (const Error& e) {
        std::cerr << "error: " << e.message << "\n";
        return Status(StatusCode::kInvalidArgument, e.message);
    }
}

// ─── Region Count ───────────────────────────────────────────────────────

std::string ParseStringArg(int argc, char** argv, const std::string& flag) {
    for (int i = 0; i < argc; ++i) {
        if (std::string(argv[i]) == flag && i + 1 < argc) return argv[i + 1];
    }
    return "";
}

uint64_t ParseU64Arg(int argc, char** argv, const std::string& flag, uint64_t def = 0) {
    auto s = ParseStringArg(argc, argv, flag);
    return s.empty() ? def : std::stoull(s);
}

/// Parse --region chr:a-b format (Contract §10.2)
/// Also falls back to --chrom + --start + --end
struct RegionSpec {
    std::string chrom;
    uint64_t a = 0;
    uint64_t b = UINT64_MAX;
    bool valid = false;
};

RegionSpec ParseRegion(int argc, char** argv) {
    RegionSpec r;
    auto region_str = ParseStringArg(argc, argv, "--region");
    if (!region_str.empty()) {
        // Parse chr:a-b
        auto colon = region_str.find(':');
        if (colon != std::string::npos) {
            r.chrom = region_str.substr(0, colon);
            auto dash = region_str.find('-', colon + 1);
            if (dash != std::string::npos) {
                r.a = std::stoull(region_str.substr(colon + 1, dash - colon - 1));
                r.b = std::stoull(region_str.substr(dash + 1));
            }
        }
        r.valid = !r.chrom.empty();
    } else {
        // Fallback: --chrom + --start + --end
        r.chrom = ParseStringArg(argc, argv, "--chrom");
        r.a = ParseU64Arg(argc, argv, "--start", 0);
        r.b = ParseU64Arg(argc, argv, "--end", UINT64_MAX);
        r.valid = !r.chrom.empty();
    }
    return r;
}

/// Parse --strand flag (overrides index strand mode)
StrandMode ParseStrandMode(int argc, char** argv, StrandMode index_default) {
    auto s = ParseStringArg(argc, argv, "--strand");
    if (s == "single" || s == "SingleStrand") return StrandMode::SingleStrand;
    if (s == "complete" || s == "StrandComplete") return StrandMode::StrandComplete;
    return index_default;
}

Status HandleRegionCount(int argc, char** argv) {
    auto index_path = ParseIndex(argc, argv);
    auto pattern = ParsePattern(argc, argv);
    auto region = ParseRegion(argc, argv);

    if (index_path.empty() || pattern.empty() || !region.valid) {
        std::cerr << "error: --index, --pattern, --region chr:a-b (or --chrom) required\n";
        return Status(StatusCode::kInvalidArgument, "Missing args");
    }

    try {
        auto idx = ReadBsi(index_path);
        StrandMode mode = ParseStrandMode(argc, argv,
                              static_cast<StrandMode>(idx.header.strand_mode));
        uint64_t count;
        if (idx.sarange.IsBuilt()) {
            // ENHANCED tier: SARange wavelet tree (Architecture §5.3)
            count = RegionalCountSARange(pattern, region.chrom, region.a, region.b,
                                          idx.fm, idx.bv.B_read, idx.bv.B_window,
                                          idx.windows, idx.reads, idx.sarange,
                                          idx.chrom_to_id, mode);
        } else if (!idx.meta_directory.empty() && !idx.map_directory.empty()) {
            count = RegionalCountLazy(pattern, region.chrom, region.a, region.b,
                                       idx.fm, idx.bv.B_read, idx.bv.B_window,
                                       idx.windows,
                                       idx.meta_payload, idx.map_payload,
                                       idx.meta_directory, idx.map_directory,
                                       idx.header.meta_codec_id, idx.header.map_codec_id,
                                       idx.chrom_to_id, mode);
        } else {
            count = RegionalCount(pattern, region.chrom, region.a, region.b,
                                   idx.fm, idx.bv.B_read, idx.bv.B_window,
                                   idx.windows, idx.reads, idx.chrom_to_id, mode);
        }
        if (HasFlag(argc, argv, "--json")) {
            std::string mode_str = (mode == StrandMode::StrandComplete) ? "strand_complete" : "single_strand";
            std::cout << "{\"region_count\":" << count
                      << ",\"mode\":\"BASE\""
                      << ",\"strand_mode\":\"" << mode_str << "\""
                      << "}\n";
        } else {
            std::cout << count << "\n";
        }
        return Status::Ok();
    } catch (const Error& e) {
        std::cerr << "error: " << e.message << "\n";
        return Status(StatusCode::kInvalidArgument, e.message);
    }
}

// ─── Region Exists ──────────────────────────────────────────────────────

Status HandleRegionExists(int argc, char** argv) {
    auto index_path = ParseIndex(argc, argv);
    auto pattern = ParsePattern(argc, argv);
    auto region = ParseRegion(argc, argv);
    uint64_t threshold = ParseU64Arg(argc, argv, "--threshold", 1);

    if (index_path.empty() || pattern.empty() || !region.valid) {
        std::cerr << "error: --index, --pattern, --region chr:a-b (or --chrom) required\n";
        return Status(StatusCode::kInvalidArgument, "Missing args");
    }

    try {
        auto idx = ReadBsi(index_path);
        StrandMode mode = ParseStrandMode(argc, argv,
                              static_cast<StrandMode>(idx.header.strand_mode));
        bool exists;
        if (!idx.meta_directory.empty() && !idx.map_directory.empty()) {
            exists = RegionalExistsLazy(pattern, threshold, region.chrom, region.a, region.b,
                                         idx.fm, idx.bv.B_read, idx.bv.B_window,
                                         idx.windows,
                                         idx.meta_payload, idx.map_payload,
                                         idx.meta_directory, idx.map_directory,
                                         idx.header.meta_codec_id, idx.header.map_codec_id,
                                         idx.chrom_to_id, mode);
        } else {
            exists = RegionalExists(pattern, threshold, region.chrom, region.a, region.b,
                                     idx.fm, idx.bv.B_read, idx.bv.B_window,
                                     idx.windows, idx.reads, idx.chrom_to_id, mode);
        }
        if (HasFlag(argc, argv, "--json")) {
            std::string mode_str = (mode == StrandMode::StrandComplete) ? "strand_complete" : "single_strand";
            std::cout << "{\"region_exists\":" << (exists ? "true" : "false")
                      << ",\"mode\":\"BASE\""
                      << ",\"strand_mode\":\"" << mode_str << "\""
                      << "}\n";
        } else {
            std::cout << (exists ? "true" : "false") << "\n";
        }
        return Status::Ok();
    } catch (const Error& e) {
        std::cerr << "error: " << e.message << "\n";
        return Status(StatusCode::kInvalidArgument, e.message);
    }
}

// ─── Reconstruct ────────────────────────────────────────────────────────

/// Parse --streams value (comma-separated: seq,qual,meta,map)
struct StreamSelection {
    bool seq = true, qual = true, meta = true, map_ = true;
};
StreamSelection ParseStreams(int argc, char** argv) {
    auto val = ParseStringArg(argc, argv, "--streams");
    if (val.empty()) return {};  // all streams
    StreamSelection sel{false, false, false, false};
    std::istringstream iss(val);
    std::string tok;
    while (std::getline(iss, tok, ',')) {
        if (tok == "seq")  sel.seq = true;
        if (tok == "qual") sel.qual = true;
        if (tok == "meta") sel.meta = true;
        if (tok == "map")  sel.map_ = true;
    }
    return sel;
}

/// Parse --read-ids as comma-separated list
std::vector<uint64_t> ParseReadIds(int argc, char** argv) {
    auto val = ParseStringArg(argc, argv, "--read-ids");
    if (val.empty()) {
        // Also check --read-id (single)
        val = ParseStringArg(argc, argv, "--read-id");
    }
    std::vector<uint64_t> ids;
    if (val.empty()) return ids;
    std::istringstream iss(val);
    std::string tok;
    while (std::getline(iss, tok, ',')) {
        ids.push_back(std::stoull(tok));
    }
    return ids;
}

Status HandleReconstruct(int argc, char** argv) {
    auto index_path = ParseIndex(argc, argv);
    auto output = ParseOutput(argc, argv, "");
    auto read_ids = ParseReadIds(argc, argv);
    auto streams = ParseStreams(argc, argv);
    bool allow_lossy = HasFlag(argc, argv, "--allow-lossy");

    if (index_path.empty()) {
        std::cerr << "error: --index <file.bsi> required\n";
        return Status(StatusCode::kInvalidArgument, "Missing --index");
    }
    if (output.empty()) {
        std::cerr << "error: --output <file.bam> required\n";
        return Status(StatusCode::kInvalidArgument, "Missing --output");
    }

    try {
        auto idx = ReadBsi(index_path);
        const auto& h = idx.header;

        // Check lossy mode
        if (h.is_lossless == 0 && !allow_lossy) {
            std::cerr << "error: Index was built with lossy encoding.\n"
                      << "       Use --allow-lossy to proceed.\n";
            return Status(StatusCode::kInvalidArgument, "Lossy index, use --allow-lossy");
        }
        if (h.is_lossless == 0) {
            std::cerr << "WARNING: Lossy index. Quality scores may differ from originals.\n";
        }

        // Determine which reads to reconstruct
        std::vector<uint64_t> target_ids;
        if (read_ids.empty()) {
            target_ids.resize(idx.reads.size());
            std::iota(target_ids.begin(), target_ids.end(), 0);
        } else {
            target_ids = read_ids;
        }

        // Validate read IDs
        for (auto rid : target_ids) {
            if (rid >= idx.reads.size()) {
                std::cerr << "error: read_id " << rid << " out of range (max: "
                          << idx.reads.size() - 1 << ")\n";
                return Status(StatusCode::kInvalidArgument, "read_id out of range");
            }
        }

        // Write BAM output via htslib
        htsFile* out_fp = hts_open(output.c_str(), "wb");
        if (!out_fp) {
            std::cerr << "error: cannot open output BAM: " << output << "\n";
            return Status(StatusCode::kInvalidArgument, "Cannot open output BAM");
        }

        // Build SAM header from chromosome table
        sam_hdr_t* out_hdr = sam_hdr_init();
        std::string hdr_text = "@HD\tVN:1.6\tSO:coordinate\n";
        for (const auto& ce : h.chrom_name_table) {
            hdr_text += "@SQ\tSN:" + ce.name + "\tLN:999999999\n";
        }
        hdr_text += "@PG\tID:bamsix\tPN:bamsix\tVN:" + std::string(h.bamsix_version) + "\n";
        sam_hdr_parse(hdr_text.size(), hdr_text.c_str());
        // Use sam_hdr_add_lines instead
        if (sam_hdr_add_lines(out_hdr, hdr_text.c_str(), hdr_text.size()) < 0) {
            std::cerr << "error: failed to build SAM header\n";
            hts_close(out_fp);
            sam_hdr_destroy(out_hdr);
            return Status(StatusCode::kInvalidArgument, "Failed to build SAM header");
        }
        if (sam_hdr_write(out_fp, out_hdr) < 0) {
            std::cerr << "error: failed to write BAM header\n";
            hts_close(out_fp);
            sam_hdr_destroy(out_hdr);
            return Status(StatusCode::kInvalidArgument, "Failed to write BAM header");
        }

        // Write each read as a BAM record
        bam1_t* rec = bam_init1();
        static const char CODE_TO_BASE[] = "ACGTN";
        uint64_t written = 0;

        for (auto rid : target_ids) {
            const auto& rd = idx.reads[rid];

            // Construct read name — use stored QNAME for lossless reconstruction
            std::string qname = rd.qname.empty()
                ? ("read_" + std::to_string(rid))  // fallback for old .bsi files
                : rd.qname;

            // Build sequence string
            std::string seq_str(rd.seq.size(), 'N');
            if (streams.seq) {
                for (size_t j = 0; j < rd.seq.size(); ++j) {
                    seq_str[j] = (rd.seq[j] < 5) ? CODE_TO_BASE[rd.seq[j]] : 'N';
                }
            }

            // Build CIGAR
            std::vector<uint32_t> cigar_bam;
            if (streams.meta) {
                cigar_bam.resize(rd.cigar.size());
                for (size_t j = 0; j < rd.cigar.size(); ++j) {
                    cigar_bam[j] = bam_cigar_gen(rd.cigar[j].len, rd.cigar[j].op);
                }
            } else {
                // Default: simple match
                cigar_bam.push_back(bam_cigar_gen(rd.seq.size(), BAM_CMATCH));
            }

            // Quality — decode from S_qual stream payload using directory
            // Contract §2.3 F4: handle both per-read and block-level directories
            std::vector<uint8_t> qual_decoded;
            if (streams.qual && !idx.qual_payload.empty()) {
                if (std::holds_alternative<StreamDirectoryBlockLevel>(idx.qual_directory)) {
                    // Block-level directory
                    const auto& block_dir = std::get<StreamDirectoryBlockLevel>(idx.qual_directory);
                    qual_decoded = DecodeQualFromBlock(idx.qual_payload, block_dir,
                                                       idx.qual_block_size, rid,
                                                       idx.header.qual_codec_id);
                } else {
                    // Per-read directory
                    const auto& per_read_dir = std::get<StreamDirectoryPerRead>(idx.qual_directory);
                    if (rid < per_read_dir.size()) {
                        qual_decoded = DecodeQualRead(idx.qual_payload, per_read_dir[rid],
                                                      idx.header.qual_codec_id);
                    }
                }
            }
            // Fallback to rd.qual if payload decode produced nothing
            // (backward compat with old .bsi files that store qual in read metadata)
            if (qual_decoded.empty() && !rd.qual.empty()) {
                qual_decoded = rd.qual;
            }
            std::string qual_str(qual_decoded.empty() ? rd.seq.size() : qual_decoded.size(), '\xFF');
            if (streams.qual) {
                for (size_t j = 0; j < qual_str.size() && j < qual_decoded.size(); ++j) {
                    qual_str[j] = static_cast<char>(qual_decoded[j]);
                }
            }

            // Set BAM fields
            bam_set1(rec,
                     qname.size(), qname.c_str(),
                     streams.meta ? rd.flag : 0,
                     streams.map_ ? static_cast<int32_t>(rd.chrom_id) : -1,
                     streams.map_ ? static_cast<hts_pos_t>(rd.pos - 1) : -1,  // 1-based → 0-based
                     rd.mapq,  // Use stored MAPQ (Contract §2.1)
                     cigar_bam.size(), cigar_bam.data(),
                     -1, -1, 0,  // mate info (not stored)
                     seq_str.size(), seq_str.c_str(),
                     qual_str.c_str(),
                     rd.aux_data.size());  // aux data length

            // Append raw aux tag bytes (Contract §2.2 lossless)
            if (streams.meta && !rd.aux_data.empty()) {
                std::memcpy(bam_get_aux(rec), rd.aux_data.data(), rd.aux_data.size());
            }

            if (sam_write1(out_fp, out_hdr, rec) < 0) {
                std::cerr << "error: failed to write BAM record for read " << rid << "\n";
                break;
            }
            ++written;
        }

        bam_destroy1(rec);
        sam_hdr_destroy(out_hdr);
        hts_close(out_fp);

        std::cerr << "[reconstruct] Wrote " << written << " records to " << output << "\n";
        return Status::Ok();
    } catch (const Error& e) {
        std::cerr << "error: " << e.message << "\n";
        return Status(StatusCode::kInvalidArgument, e.message);
    }
}

}  // namespace

Status DispatchKnownCommand(const std::string& cmd) {
    if (cmd == "build") {
        std::cerr << "Usage: bamsix build <input.bam> [-o <output.bsi>]\n";
        return Status::Ok();
    }

    static const std::unordered_set<std::string> kKnownCommands = {
        "build", "count", "exists", "locate",
        "region-count", "region-exists",
        "reconstruct", "info", "verify",
    };

    if (kKnownCommands.contains(cmd)) {
        std::cerr << "error: subcommand not implemented yet: " << cmd << "\n";
        return Status(StatusCode::kInvalidArgument,
                      "subcommand not implemented yet: " + cmd);
    }

    std::cerr << "error: unknown subcommand: " << cmd << "\n";
    std::cerr << "usage: bamsix <subcommand>\n";
    std::cerr << "subcommands: version, build, count, exists, locate, "
              << "region-count, region-exists, reconstruct, info, verify\n";
    return Status(StatusCode::kInvalidArgument, "unknown subcommand: " + cmd);
}

Status DispatchCommand(const std::string& cmd, int argc, char** argv) {
    // ── Contract §4.5 / §0.8: Refuse approximate-match queries with
    //    structured NOT_IMPLEMENTED_V1 error. v1.0 is exact-match only.
    static const std::unordered_set<std::string> kQueryCmds = {
        "count", "exists", "locate", "region-count", "region-exists",
    };
    static const std::unordered_set<std::string> kApproxFlags = {
        "--mismatch", "--mismatches", "--k-mismatch", "--k-edit",
        "--edit-distance", "--approximate", "--approx", "--hamming",
    };
    if (kQueryCmds.contains(cmd)) {
        for (int i = 2; i < argc; ++i) {
            if (kApproxFlags.contains(std::string(argv[i]))) {
                std::cerr << "error: Approximate matching (--" << argv[i]
                          << ") is not implemented in BAMSIX v1.0.\n"
                          << "       v1.0 supports exact-match queries only.\n"
                          << "       Approximate matching is planned for v2.0.\n";
                return Status(StatusCode::kInvalidArgument,
                              "NOT_IMPLEMENTED_V1: approximate matching not available");
            }
        }
    }

    if (cmd == "build")         return HandleBuild(argc, argv);
    if (cmd == "count")         return HandleCount(argc, argv);
    if (cmd == "exists")        return HandleExists(argc, argv);
    if (cmd == "locate")        return HandleLocate(argc, argv);
    if (cmd == "verify")        return HandleVerify(argc, argv);
    if (cmd == "info")          return HandleInfo(argc, argv);
    if (cmd == "region-count")  return HandleRegionCount(argc, argv);
    if (cmd == "region-exists") return HandleRegionExists(argc, argv);
    if (cmd == "reconstruct")   return HandleReconstruct(argc, argv);

    return DispatchKnownCommand(cmd);
}

}  // namespace bamsix
