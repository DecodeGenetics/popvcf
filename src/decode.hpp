#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include <popvcf/constants.hpp>

namespace popvcf
{
//! Buffer size when decoding
long constexpr DEC_BUFFER_SIZE{65536};

//! Data type of an array buffer
using Tdec_array_buf = std::array<char, DEC_BUFFER_SIZE>;

class DecodeData
{
public:
  std::size_t bytes_read{0};
  std::size_t remaining_bytes{0};
  std::size_t field{0}; // current vcf field
  std::size_t b{0};     // begin index in buffer_in
  std::size_t i{b};     // index in buffer_in
  std::size_t o{0};     // output index

  std::vector<std::string> unique_fields{};
};

//! Decodes an input buffer. Output is written in \a buffer_out .
template <typename Tbuffer_out>
void decode_buffer(Tbuffer_out & buffer_out, Tdec_array_buf & buffer_in, DecodeData & ed);

//! Decode an encoded popVCF
void decode_file(std::string const & popvcf_fn, bool const is_bgzf_input);

} // namespace popvcf
