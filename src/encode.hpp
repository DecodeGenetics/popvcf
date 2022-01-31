#pragma once

#include <cstdint>
#include <string>

#include <parallel_hashmap/phmap.h>

namespace popvcf
{
//! Buffer size when encoding
long constexpr ENC_BUFFER_SIZE{4 * 65536};

//! Data type of an array buffer
using Tarray_buf = std::array<char, ENC_BUFFER_SIZE>; //!< Buffer type

class EncodeData
{
public:
  std::size_t bytes_read{0};
  std::size_t remaining_bytes{0};
  std::size_t field{0};   // current vcf field
  std::size_t b{0};       // begin index in buffer_in
  std::size_t i{b};       // index in buffer_in
  std::size_t o{0};       // output index
  bool header_line{true}; //!< True iff in header line

  uint32_t num_unique_fields{0};
  phmap::flat_hash_map<std::string, uint32_t> map_to_unique_fields{};

  inline void clear_line()
  {
    field = 0;                    // reset field index
    map_to_unique_fields.clear(); // clear map every line
    num_unique_fields = 0;        // clear the number of unique fields
  }
};

//! Encodes an input buffer. Output is written in \a buffer_out .
template <typename Tbuffer_out>
void encode_buffer(Tbuffer_out & buffer_out, Tarray_buf & buffer_in, EncodeData & ed);

//! Encode a gzipped file and write to stdout
void encode_file(std::string const & input_fn,
                 bool const is_bgzf_input,
                 std::string const & output_fn,
                 std::string const & output_mode,
                 bool const is_bgzf_output,
                 int const compression_threads);

} // namespace popvcf
