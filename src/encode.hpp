#pragma once

#include <cstdint>
#include <string>
#include <vector>

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

  std::string prev_contig{};
  uint64_t prev_pos{0};
  std::vector<std::string> prev_unique_fields{};
  //std::vector<uint32_t> prev_field2uid{};
  phmap::flat_hash_map<std::string, uint32_t> prev_map_to_unique_fields{};

  std::string contig{};
  uint64_t pos{0};
  std::vector<std::string> unique_fields{};
  //std::vector<uint32_t> field2uid{};
  phmap::flat_hash_map<std::string, uint32_t> map_to_unique_fields{};

  inline void clear_line()
  {
    field = 0;                    // reset field index

    std::swap(prev_contig, contig);
    prev_pos = pos;
    std::swap(prev_unique_fields, unique_fields);
    //std::swap(prev_field2uid, field2uid);
    std::swap(prev_map_to_unique_fields, map_to_unique_fields);

    contig.resize(0);
    pos = 0;
    unique_fields.resize(0);
    //field2uid.resize(0);
    map_to_unique_fields.clear(); // clear map every line
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
