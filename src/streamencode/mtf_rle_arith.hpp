#pragma once
/// BWT → MTF → RLE → Arithmetic codec (Contract §2.4)
#include <cstdint>
#include <vector>

namespace bamsix {

/// Move-to-Front transform over alphabet size sigma.
std::vector<uint8_t> MtfEncode(const std::vector<uint8_t>& input, uint8_t sigma);
std::vector<uint8_t> MtfDecode(const std::vector<uint8_t>& input, uint8_t sigma);

/// Run-Length Encoding of MTF output.
/// Format: for each run, (symbol:1, length:varint).
std::vector<uint8_t> RleEncode(const std::vector<uint8_t>& input);
std::vector<uint8_t> RleDecode(const std::vector<uint8_t>& encoded, uint64_t expected_len);

/// 0th-order adaptive arithmetic coder.
std::vector<uint8_t> ArithEncode(const std::vector<uint8_t>& input);
std::vector<uint8_t> ArithDecode(const std::vector<uint8_t>& encoded, uint64_t expected_len);

}  // namespace bamsix
