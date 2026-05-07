/*
 * gen_test_bam.c — Generate a 10-read synthetic BAM for V1 gate test.
 * Uses htslib's bam_set1() for correct record construction.
 */
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    const char* name;
    const char* seq;
    const char* qual;
    int pos;          /* 0-based */
    const char* cigar_str;
    uint16_t flag;
} TestRead;

/* 10 reads with CIGARs whose query-consuming ops sum to |seq| */
static TestRead reads[10] = {
    {"read01", "ACGTACGTACGTACGT", "IIIIIIIIIIIIIIII", 100, "16M", 0},
    {"read02", "CCGGTTAACCGGTTAA", "HHHHHHHHHHHHHHHH", 200, "16M", 0},
    {"read03", "AAAACCCCGGGGTTTT", "EEEEEEEEEEEEEEEE", 300, "16M", 0},
    {"read04", "ACGTACGTACGT",     "IIIIIIIIIIII",     400, "4M2I6M",  0},   /* 4+2+6=12 ✓ */
    {"read05", "TTTTCCCCAAAA",     "HHHHHHHHHHHH",     500, "6M2D6M",  0},   /* 6+6=12 ✓ */
    {"read06", "GGGGCCCCAAAA",     "EEEEEEEEEEEE",     600, "2S8M2S",  0},   /* 2+8+2=12 ✓ */
    {"read07", "ACGTACGTACGTACGT", "IIIIIIIIIIIIIIII", 700, "16M", 16},      /* reverse */
    {"read08", "AACCGGTTAACCGGTT", "HHHHHHHHHHHHHHHH", 800, "8M1I7M",  0},   /* 8+1+7=16 ✓ */
    {"read09", "CCCCAAAATTTTGGGG", "EEEEEEEEEEEEEEEE", 900, "3S10M3S", 0},   /* 3+10+3=16 ✓ */
    {"read10", "ACGTACGT",         "IIIIIIII",         1000, "8M",     0},
};

/* Parse CIGAR string to uint32_t array */
static int parse_cigar(const char* s, uint32_t* cigar, int max_ops) {
    int n = 0;
    const char* p = s;
    while (*p && n < max_ops) {
        int len = 0;
        while (*p >= '0' && *p <= '9') { len = len * 10 + (*p - '0'); p++; }
        int op;
        switch (*p) {
            case 'M': op = BAM_CMATCH; break;
            case 'I': op = BAM_CINS; break;
            case 'D': op = BAM_CDEL; break;
            case 'N': op = BAM_CREF_SKIP; break;
            case 'S': op = BAM_CSOFT_CLIP; break;
            case 'H': op = BAM_CHARD_CLIP; break;
            case '=': op = BAM_CEQUAL; break;
            case 'X': op = BAM_CDIFF; break;
            default: return -1;
        }
        cigar[n++] = bam_cigar_gen(len, op);
        p++;
    }
    return n;
}

int main(int argc, char** argv) {
    const char* outpath = "data/test/synthetic_10reads.bam";
    if (argc > 1) outpath = argv[1];

    htsFile* fp = hts_open(outpath, "wb");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", outpath); return 1; }

    sam_hdr_t* hdr = sam_hdr_init();
    sam_hdr_add_line(hdr, "HD", "VN", "1.6", "SO", "coordinate", NULL);
    sam_hdr_add_line(hdr, "SQ", "SN", "chr1", "LN", "10000", NULL);
    if (sam_hdr_write(fp, hdr) < 0) { fprintf(stderr, "Header write fail\n"); return 1; }

    bam1_t* b = bam_init1();

    for (int i = 0; i < 10; i++) {
        TestRead* r = &reads[i];
        int seq_len = (int)strlen(r->seq);

        uint32_t cigar[16];
        int n_cigar = parse_cigar(r->cigar_str, cigar, 16);
        if (n_cigar < 0) { fprintf(stderr, "Bad CIGAR for read %d\n", i); return 1; }

        /* Use bam_set1 for correct record construction */
        int ret = bam_set1(b,
            strlen(r->name),        /* l_qname (excluding NUL) */
            r->name,                /* qname */
            r->flag,                /* flag */
            0,                      /* tid (chr1 = 0) */
            r->pos,                 /* pos (0-based) */
            60,                     /* mapq */
            n_cigar,                /* n_cigar */
            cigar,                  /* cigar */
            -1,                     /* mtid */
            -1,                     /* mpos */
            0,                      /* isize */
            seq_len,                /* l_seq */
            r->seq,                 /* seq (ASCII) */
            r->qual,               /* qual (ASCII Phred+33) */
            0                       /* l_aux */
        );
        if (ret < 0) { fprintf(stderr, "bam_set1 failed for read %d\n", i); return 1; }

        if (sam_write1(fp, hdr, b) < 0) {
            fprintf(stderr, "Write failed for read %d\n", i);
            return 1;
        }
    }

    bam_destroy1(b);
    sam_hdr_destroy(hdr);
    hts_close(fp);
    printf("Wrote %s with 10 reads\n", outpath);
    return 0;
}
