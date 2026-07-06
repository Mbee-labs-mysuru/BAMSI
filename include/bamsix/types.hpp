#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <variant>
#include <vector>

namespace bamsix {

// ─── Alphabet Encoding (Architecture §0.2 / §2.3) ──────────────────────────

constexpr uint8_t CODE_A     = 0;
constexpr uint8_t CODE_C     = 1;
constexpr uint8_t CODE_G     = 2;
constexpr uint8_t CODE_T     = 3;
constexpr uint8_t CODE_N     = 4;
constexpr uint8_t CODE_SEP   = 5;  // '#' separator
constexpr uint8_t CODE_SENT  = 6;  // '$' sentinel — conceptual only, never stored
constexpr uint8_t SIGMA      = 6;  // alphabet size for stored codes {0..5}
constexpr uint8_t SIGMA_EXT  = 7;  // extended alphabet including sentinel

// ─── CIGAR operations ────────────────────────────────────────────────────────

struct CigarOp {
    uint8_t  op;   // BAM CIGAR op code: 0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 6=P, 7==, 8=X
    uint32_t len;
};

using CigarRecord = std::vector<CigarOp>;

// ─── Ingestion (Architecture §4.1) ──────────────────────────────────────────

struct RawRead {
    std::vector<uint8_t> seq;          // encoded bases: codes 0..4
    std::string          chrom;        // chromosome name string
    std::string          qname;        // read name (QNAME) — required for lossless reconstruct
    uint32_t             flag;         // BAM FLAG field
    uint64_t             pos;          // 1-based SAM POS (converted once at ingestion)
    uint8_t              mapq;         // mapping quality — required for lossless reconstruct
    CigarRecord          cigar;        // ordered list of (op, len) pairs
    std::vector<uint8_t> qual;         // Phred scores in [0,93]; |qual| == |seq|
    uint32_t             source_file_id;   // 0-based index in input BAM list
    uint64_t             bam_offset;       // stable record index within source file
    std::vector<uint8_t> aux_data;         // raw BAM aux tag bytes (Contract §2.2)
};

// ─── Ordering (Architecture §4.2) ───────────────────────────────────────────

struct OrderedRead : RawRead {
    uint32_t chrom_id;    // 0-based lexicographic rank of chrom name
    uint64_t read_id;     // 0-based rank under total order ≺
};

// ─── Sequence Builder (Architecture §4.3) ────────────────────────────────────

struct SequenceBundle {
    std::vector<uint8_t>  S;            // r0#r1#...#r_{N-1}, codes 0..5
    std::vector<uint64_t> readStarts;   // readStarts[i] = start of read i in S
};

// ─── Stream Directories (Architecture §3 / §7.6) ────────────────────────────

struct StreamDirectoryEntry {
    uint64_t offset;   // byte offset of this read's block in the stream payload
    uint32_t length;   // byte length of this read's block
};

using StreamDirectoryPerRead = std::vector<StreamDirectoryEntry>;

struct BlockDirectoryEntry {
    uint64_t block_offset;
    uint32_t block_length;
    uint32_t first_read_id;
};
using StreamDirectoryBlockLevel = std::vector<BlockDirectoryEntry>;

using StreamDirectory = std::variant<StreamDirectoryPerRead, StreamDirectoryBlockLevel>;

struct Directories {
    // MANDATORY per-read directories (hot-path for query-time mapping)
    StreamDirectoryPerRead  meta;
    StreamDirectoryPerRead  map;

    // OPTIONAL granularity (per-read OR block-level)
    StreamDirectory         seq;
    StreamDirectory         qual;
    uint32_t                seq_block_size;   // 0 = per-read
    uint32_t                qual_block_size;  // 0 = per-read
};

// ─── Encoded Streams ─────────────────────────────────────────────────────────

struct EncodedStreams {
    std::vector<uint8_t> S_seq;
    std::vector<uint8_t> S_qual;
    std::vector<uint8_t> S_meta;
    std::vector<uint8_t> S_map;
    Directories          dirs;
};

// ─── FM-Index (Architecture §4.6) ────────────────────────────────────────────

struct FMIndex {
    std::vector<uint8_t>   BWT;          // stored BWT payload, non-sentinel rows
    std::array<uint64_t, 7> C;           // C[a] = # SA rows with suffix < a
    // Occ structure is stored separately (wavelet tree / RRR)
    uint64_t               sample_step_s;
    std::vector<uint64_t>  SA_samples;
    uint64_t               sentinel_row;

    // Optional ISA samples
    bool                   has_isa_samples = false;
    uint64_t               sample_step_s_prime = 0;
    std::vector<uint64_t>  ISA_samples;
};

// ─── SA interval from backward search ────────────────────────────────────────

struct SAInterval {
    uint64_t lo;  // inclusive
    uint64_t hi;  // exclusive
    bool empty() const { return lo >= hi; }
    uint64_t size() const { return (lo < hi) ? (hi - lo) : 0; }
};

// ─── Windows (Architecture §4.8) ─────────────────────────────────────────────

struct Window {
    uint32_t chrom_id;
    uint64_t l;                  // inclusive start in S
    uint64_t r;                  // inclusive end in S
    uint64_t first_read_id;
    uint64_t last_read_id;
    uint64_t genomic_start;      // 1-based
    uint64_t genomic_end;        // 1-based
};

using WindowTable = std::vector<Window>;

// ─── Query Types (Architecture §6) ───────────────────────────────────────────

enum class QueryStrand : uint8_t { Forward, Reverse };

struct MappingResult {
    uint32_t    chrom_id;
    uint64_t    p_min;        // 1-based inclusive
    uint64_t    p_max;        // 1-based inclusive
    uint64_t    read_id;
    QueryStrand query_strand;
};

struct Match {
    std::string chrom;
    uint64_t    p_min;
    uint64_t    p_max;
    uint64_t    read_id;
    QueryStrand query_strand;
    uint64_t    sa_row = 0;   // Contract §4.1: 6th tie-break key for sorted mode
};

// ─── Build Configuration (Architecture §3) ───────────────────────────────────

enum class StrandMode : uint8_t { StrandComplete = 0, SingleStrand = 1 };

/// OverlapMode controls how overlapping pattern occurrences are counted.
/// All  — every occurrence position in S (FM-index default; overlapping matches)
/// None — greedy non-overlapping: after a match at position p, skip to p+|P|
///         (grep -o compatible semantics)
enum class OverlapMode : uint8_t { All = 0, None = 1 };

enum class QualCodec : uint8_t {
    RANGE_CYCLE   = 0x01,
    RANS_DELTA    = 0x02,
    ZSTD_DICT     = 0x03,
    BINNED_RANGE  = 0x04,
};

enum class MetaCodec : uint8_t {
    TYPED_SPLIT   = 0x01,
    ZSTD_FALLBACK = 0x02,
};

enum class MapCodec : uint8_t {
    DELTA_RANGE   = 0x01,
    RAW           = 0x02,
};

struct BuildConfig {
    uint64_t   window_size_T         = 100000;
    uint64_t   sample_step_s         = 64;
    uint64_t   sample_step_s_prime   = 0;       // 0 = ISA disabled
    uint8_t    entropy_order_k       = 6;
    StrandMode strand_mode           = StrandMode::StrandComplete;
    bool       enable_sarange        = false;
    bool       shared_bwt            = true;

    QualCodec  qual_codec            = QualCodec::RANGE_CYCLE;
    uint8_t    qual_lossy_bins       = 0;       // 0 = lossless
    MetaCodec  meta_codec            = MetaCodec::TYPED_SPLIT;
    MapCodec   map_codec             = MapCodec::DELTA_RANGE;

    uint32_t   seq_block_size        = 1024;
    uint32_t   qual_block_size       = 1024;

    bool       enable_bidirectional  = false;
    uint8_t    recommended_seed_length = 16;

    bool       allow_parallel_sa     = false;
    uint32_t   build_threads         = 0;

    std::string reference_fasta_path;
    std::array<uint8_t, 32> reference_sha256 = {};
};

// ─── Structured Errors (Architecture §8.1) ───────────────────────────────────

enum class ErrorCode : uint16_t {
    OK = 0,
    INVALID_BAM_INPUT,
    CORRUPT_BSI,
    UNSUPPORTED_CIGAR_OP,
    INVALID_PATTERN,
    EMPTY_PATTERN,
    SEPARATOR_POSITION,
    CHECKSUM_MISMATCH,
    ORDERING_HASH_MISMATCH,
    MANIFEST_MISMATCH,
    VERSION_MISMATCH,
    BUILD_VALIDATION_FAILED,
    STREAM_DECODE_ERROR,
    NOT_IMPLEMENTED_V1,
    REFERENCE_MISMATCH,
    LOSSY_RECONSTRUCTION,
    UNSUPPORTED_CODEC,
};

struct Error {
    ErrorCode   code;
    std::string message;
};

// ─── BSI file header (Architecture §7.1) ─────────────────────────────────────

struct ChromEntry {
    uint32_t chrom_id;
    std::string name;
};

struct BsiHeader {
    // Magic
    char     magic[4]       = {'B','S','I','X'};
    uint16_t version        = 6;

    // Build provenance
    char     bamsix_version[16] = {};
    uint8_t  host_os_id     = 0;
    uint8_t  cpu_arch_id    = 0;
    uint64_t build_timestamp_utc = 0;
    uint8_t  is_lossless    = 1;

    // Input manifest
    uint32_t source_file_count = 0;
    std::array<uint8_t, 32> source_manifest_hash = {};
    std::array<uint8_t, 32> ordering_hash = {};

    // Dimensions
    uint64_t S_length       = 0;
    uint64_t N_reads        = 0;
    uint32_t N_windows      = 0;

    // Sampling
    uint32_t sample_step_s  = 64;
    uint8_t  has_isa_samples = 0;
    uint32_t sample_step_s_prime = 0;

    // Tier
    uint8_t  enable_sarange = 0;
    uint8_t  sarange_variant = 0;
    uint8_t  shared_bwt     = 1;
    uint8_t  enable_bidirectional = 0;
    uint8_t  recommended_seed_length = 16;

    // Window / entropy
    uint64_t window_size_T  = 100000;
    uint8_t  entropy_order_k = 6;

    // Codecs
    uint8_t  qual_codec_id  = 0x01;
    uint8_t  qual_lossy_bins = 0;
    uint8_t  meta_codec_id  = 0x01;
    uint8_t  map_codec_id   = 0x01;

    // Strand / sentinel
    uint8_t  strand_mode    = 0;
    uint64_t sentinel_row   = 0;

    // Chromosome table
    uint32_t chrom_count    = 0;
    std::vector<ChromEntry> chrom_name_table;

    // Block sizes
    uint32_t seq_block_size  = 1024;
    uint32_t qual_block_size = 1024;

    // Parallelism / reference
    uint8_t  allow_parallel_sa = 0;
    uint8_t  reference_based_encoding = 0;
    std::array<uint8_t, 32> reference_sha256 = {};
    uint32_t flags          = 0;
};

}  // namespace bamsix
