# BAMSI Benchmark Results — v1.0.0

Generated from 20 data points.

## Compression Results

| Dataset | Tool | Build Time (s) | Compressed Size (bytes) |
|---------|------|---------------|------------------------|
| synthetic_10 | bamsi | 0.007 ± 0.001 | 2331 |

## Query Latency Results

| Dataset | Tool | Pattern | Latency (s) |
|---------|------|---------|-------------|
| synthetic_10 | bamsi | AAAA | 0.005411 ± 0.000090 |
| synthetic_10 | bamsi | ACGT | 0.004954 ± 0.000173 |
| synthetic_10 | bamsi | ACGTACGT | 0.005542 ± 0.001188 |
| synthetic_10 | bamsi | CCCC | 0.004892 ± 0.000973 |
| synthetic_10 | bamsi | GGGG | 0.005333 ± 0.001037 |
| synthetic_10 | bamsi | TGCATGCA | 0.007110 ± 0.003119 |
| synthetic_10 | bamsi | TTTT | 0.005093 ± 0.000969 |

## Verify & Reconstruct Results

| Dataset | Tool | Operation | Time (s) |
|---------|------|-----------|----------|
| synthetic_10 | bamsi | reconstruct | 0.005 ± 0.000 |
| synthetic_10 | bamsi | verify | 0.005 ± 0.000 |
