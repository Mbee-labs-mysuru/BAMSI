#include "bamsi/cli/dispatch.hpp"
#include "build.hpp"
#include "../query/query.hpp"
#include "../fmindex/fmindex.hpp"
#include "../bitvectors/bitvectors.hpp"
#include "../format/format.hpp"
#include "bamsi/config.hpp"

#include <htslib/sam.h>
#include <htslib/hts.h>

#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>
#include <sstream>
#include <iomanip>
#include <numeric>

namespace bamsi {

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
        "--read-ids", "--read-id",
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

    if (index_path.empty()) {
        std::cerr << "error: --index <file.bsi> required\n";
        return Status(StatusCode::kInvalidArgument, "Missing --index");
    }
    if (pattern.empty()) {
        std::cerr << "error: --pattern <ACGT...> required\n";
        return Status(StatusCode::kInvalidArgument, "Missing --pattern");
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
        uint64_t count = GlobalCount(pattern, idx.fm, mode);

        if (HasFlag(argc, argv, "--json")) {
            std::cout << "{\"count\":" << count << "}\n";
        } else {
            std::cout << count << "\n";
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

        auto matches = Locate(pattern, idx.fm, idx.bv.B_read,
                              idx.reads, idx.chrom_names, mode, sorted);

        if (HasFlag(argc, argv, "--json")) {
            std::cout << "{\"matches\":[";
            for (size_t i = 0; i < matches.size(); ++i) {
                if (i > 0) std::cout << ",";
                const auto& m = matches[i];
                std::cout << "{\"strand\":\"" << (m.query_strand == QueryStrand::Forward ? "+" : "-")
                          << "\",\"chrom\":\"" << m.chrom
                          << "\",\"p_min\":" << m.p_min
                          << ",\"p_max\":" << m.p_max
                          << ",\"read_id\":" << m.read_id << "}";
            }
            std::cout << "],\"count\":" << matches.size() << "}\n";
        } else {
            // TSV header
            std::cout << "strand\tchrom\tp_min\tp_max\tread_id\n";
            for (const auto& m : matches) {
                std::cout << (m.query_strand == QueryStrand::Forward ? "+" : "-")
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

    if (index_path.empty()) {
        std::cerr << "error: --index <file.bsi> required\n";
        return Status(StatusCode::kInvalidArgument, "Missing --index");
    }

    try {
        bool ok = VerifyBsi(index_path);
        if (HasFlag(argc, argv, "--json")) {
            std::cout << "{\"valid\":" << (ok ? "true" : "false") << "}\n";
        } else {
            if (ok) {
                std::cout << "PASS: " << index_path << " integrity verified\n";
            } else {
                std::cout << "FAIL: " << index_path << " checksum mismatch\n";
            }
        }
        return ok ? Status::Ok()
                  : Status(StatusCode::kInvalidArgument, "Checksum mismatch");
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
            std::cout << "  \"bamsi_version\": \"" << h.bamsi_version << "\",\n";
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
            std::cout << "]\n";
            std::cout << "}\n";
        } else {
            std::cout << "BAMSI Index: " << index_path << "\n";
            std::cout << "  Format version:      " << h.version << "\n";
            std::cout << "  BAMSI version:       " << h.bamsi_version << "\n";
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
        uint64_t count = RegionalCount(pattern, region.chrom, region.a, region.b,
                                        idx.fm, idx.bv.B_read, idx.bv.B_window,
                                        idx.windows, idx.reads, idx.chrom_to_id, mode);
        if (HasFlag(argc, argv, "--json")) {
            std::cout << "{\"region_count\":" << count << "}\n";
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
        bool exists = RegionalExists(pattern, threshold, region.chrom, region.a, region.b,
                                      idx.fm, idx.bv.B_read, idx.bv.B_window,
                                      idx.windows, idx.reads, idx.chrom_to_id, mode);
        if (HasFlag(argc, argv, "--json")) {
            std::cout << "{\"region_exists\":" << (exists ? "true" : "false") << "}\n";
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
        hdr_text += "@PG\tID:bamsi\tPN:bamsi\tVN:" + std::string(h.bamsi_version) + "\n";
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

            // Construct read name
            std::string qname = "read_" + std::to_string(rid);

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

            // Quality
            std::string qual_str(rd.qual.size(), '\xFF');
            if (streams.qual) {
                for (size_t j = 0; j < rd.qual.size(); ++j) {
                    qual_str[j] = static_cast<char>(rd.qual[j]);
                }
            }

            // Set BAM fields
            bam_set1(rec,
                     qname.size(), qname.c_str(),
                     streams.meta ? rd.flag : 0,
                     streams.map_ ? static_cast<int32_t>(rd.chrom_id) : -1,
                     streams.map_ ? static_cast<hts_pos_t>(rd.pos - 1) : -1,  // 1-based → 0-based
                     255,  // MAPQ
                     cigar_bam.size(), cigar_bam.data(),
                     -1, -1, 0,  // mate info (not stored)
                     seq_str.size(), seq_str.c_str(),
                     qual_str.c_str(),
                     0);  // no aux data

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
        std::cerr << "Usage: bamsi build <input.bam> [-o <output.bsi>]\n";
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
    std::cerr << "usage: bamsi <subcommand>\n";
    std::cerr << "subcommands: version, build, count, exists, locate, "
              << "region-count, region-exists, reconstruct, info, verify\n";
    return Status(StatusCode::kInvalidArgument, "unknown subcommand: " + cmd);
}

Status DispatchCommand(const std::string& cmd, int argc, char** argv) {
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

}  // namespace bamsi
