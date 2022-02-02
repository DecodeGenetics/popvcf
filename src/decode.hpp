#pragma once

#include <array>
#include <cassert>
#include <charconv>
#include <cstdint>
#include <string>
#include <vector>

#include "sequence_utils.hpp"

#include <popvcf/constants.hpp>

namespace popvcf
{
//! Buffer size when decoding
long constexpr DEC_BUFFER_SIZE{8 * 65536};

//! Data type of an array buffer
using Tdec_array_buf = std::array<char, DEC_BUFFER_SIZE>;

class DecodeData
{
public:
  std::size_t bytes_read{0};
  std::size_t remaining_bytes{0};
  std::size_t field{0};   // current vcf field
  std::size_t b{0};       // begin index in buffer_in
  std::size_t i{b};       // index in buffer_in
  std::size_t o{0};       // output index
  bool header_line{true}; //!< True iff in header line

  std::vector<std::string> prev_unique_fields{};
  std::vector<std::string> unique_fields{};

  inline void clear_line()
  {
    field = 0;
    std::swap(prev_unique_fields, unique_fields);
    unique_fields.resize(0);
  }
};

//! Decodes an input buffer. Output is written in \a buffer_out .
template <typename Tbuffer_out, typename Tbuffer_in>
inline void decode_buffer(Tbuffer_out & buffer_out, Tbuffer_in & buffer_in, DecodeData & dd)
{
  std::size_t constexpr N_FIELDS_SITE_DATA{9};

  // inner loop - Loops over each character in the input buffer
  while (dd.i < (dd.bytes_read + dd.remaining_bytes))
  {
    char const b_in = buffer_in[dd.i];

    if (b_in != '\t' && b_in != '\n')
    {
      ++dd.i; // we are in a vcf field
      continue;
    }

    if (dd.field == 0)
      dd.header_line = buffer_in[dd.b] == '#'; // check if in header line

    if (dd.header_line || dd.field < N_FIELDS_SITE_DATA)
    {
      // write field without any encoding
      ++dd.i; // adds '\t' or '\n'
      std::copy(&buffer_in[dd.b], &buffer_in[dd.i], std::back_inserter(buffer_out));
    }
    else
    {
      if (buffer_in[dd.b] == '%')
      {
        // unique field within the line but was seen in the previous line
        ++dd.b; // Get over '%'
        uint32_t const prev_unique_index = ascii_cstring_to_int(&buffer_in[dd.b], &buffer_in[dd.i++]);
        assert(prev_unique_index < dd.prev_unique_fields.size());
        std::string const & prior_field = dd.prev_unique_fields[prev_unique_index];
        std::copy(prior_field.begin(), prior_field.end(), std::back_inserter(buffer_out));
        dd.unique_fields.emplace_back(prior_field);
        buffer_out.push_back(b_in);
      }
      else if (buffer_in[dd.b] >= ':')
      {
        uint32_t const unique_index = ascii_cstring_to_int(&buffer_in[dd.b], &buffer_in[dd.i++]);
        assert(unique_index < dd.unique_fields.size());
        std::string const & prior_field = dd.unique_fields[unique_index];
        std::copy(prior_field.begin(), prior_field.end(), std::back_inserter(buffer_out));
        buffer_out.push_back(b_in);
      }
      else
      {
        // add a new unique field and write field without any encoding
        dd.unique_fields.emplace_back(&buffer_in[dd.b], dd.i - dd.b);
        std::copy(buffer_in.begin() + dd.b, buffer_in.begin() + (++dd.i), std::back_inserter(buffer_out));
      }
    }

    assert(b_in == buffer_in[dd.i - 1]);
    dd.b = dd.i;

    if (b_in == '\n')
      dd.clear_line();
    else
      ++dd.field;
  } // ends inner loop

  // write data to the beginning of the input buffer
  dd.remaining_bytes = dd.i - dd.b;
  std::copy(&buffer_in[dd.b], &buffer_in[dd.b + dd.remaining_bytes], &buffer_in[0]);
  dd.b = 0;
  dd.i = dd.remaining_bytes;
}

//! Decode an encoded popVCF
void decode_file(std::string const & popvcf_fn, bool const is_bgzf_input);

} // namespace popvcf
