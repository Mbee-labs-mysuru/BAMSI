/// Generate a 10,000-read synthetic BAM for V2 gate testing.
/// Reads are random sequences of length 50-150bp, mapped to chr1-chr5.

#include <htslib/sam.h>
#include <htslib/hts.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdio.h>

#define N_READS 10000
#define N_CHROMS 5
#define MIN_LEN 50
#define MAX_LEN 150

static const char* chroms[N_CHROMS] = {"chr1", "chr2", "chr3", "chr4", "chr5"};
static const int chrom_lens[N_CHROMS] = {248956422, 242193529, 198295559, 190214555, 181538259};

int main(int argc, char** argv) {
    const char* outpath = (argc > 1) ? argv[1] : "data/test/synthetic_10k.bam";

    srand(42);  // deterministic

    htsFile* out = hts_open(outpath, "wb");
    if (!out) { fprintf(stderr, "Cannot open %s\n", outpath); return 1; }

    sam_hdr_t* hdr = sam_hdr_init();
    for (int i = 0; i < N_CHROMS; i++) {
        char line[256];
        snprintf(line, sizeof(line), "@SQ\tSN:%s\tLN:%d\n", chroms[i], chrom_lens[i]);
        sam_hdr_add_lines(hdr, line, 0);
    }
    if (sam_hdr_write(out, hdr) < 0) {
        fprintf(stderr, "Failed to write header\n"); return 1;
    }

    const char bases[] = "ACGT";

    for (int i = 0; i < N_READS; i++) {
        int read_len = MIN_LEN + (rand() % (MAX_LEN - MIN_LEN + 1));
        int chrom_idx = rand() % N_CHROMS;
        int pos = rand() % (chrom_lens[chrom_idx] - read_len - 1);

        // Generate sequence
        char* seq = (char*)malloc(read_len + 1);
        for (int j = 0; j < read_len; j++) {
            seq[j] = bases[rand() % 4];
        }
        seq[read_len] = '\0';

        // Quality scores (all Q30)
        char* qual = (char*)malloc(read_len + 1);
        memset(qual, 30 + 33, read_len);
        qual[read_len] = '\0';

        // CIGAR: simple match
        uint32_t cigar = bam_cigar_gen(read_len, BAM_CMATCH);

        char qname[32];
        snprintf(qname, sizeof(qname), "read_%05d", i);

        bam1_t* b = bam_init1();
        if (bam_set1(b,
                     strlen(qname), qname,
                     0,                    // FLAG = mapped, forward
                     chrom_idx, pos,
                     60,                   // MAPQ
                     1, &cigar,            // CIGAR
                     -1, -1, 0,            // mate info
                     read_len, seq, qual,
                     0) < 0) {
            fprintf(stderr, "bam_set1 failed for read %d\n", i);
            bam_destroy1(b);
            free(seq); free(qual);
            continue;
        }

        if (sam_write1(out, hdr, b) < 0) {
            fprintf(stderr, "sam_write1 failed for read %d\n", i);
        }

        bam_destroy1(b);
        free(seq);
        free(qual);
    }

    sam_hdr_destroy(hdr);
    hts_close(out);
    printf("Generated %d reads → %s\n", N_READS, outpath);
    return 0;
}
