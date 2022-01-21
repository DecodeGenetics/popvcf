#include "decode.hpp"

#include <algorithm> //std::copy
#include <array>     // std::array
#include <cstring>   // std::memmove
#include <iostream>  // std::cerr
#include <string>    // std::string
#include <vector>    // std::vector

#include "sequence_utils.hpp" // int_to_ascii

namespace
{
long constexpr DEC_BUFFER_SIZE{65536};                 //!< Buffer input size when decoding
using Tdec_buf_in = std::array<char, DEC_BUFFER_SIZE>; //!< Buffer in type
using Tdec_buf_out = std::vector<char>;                //!< Buffer out type
} // namespace

namespace popvcf
{
void decode()
{
  Tdec_buf_in buffer_in;   // input buffer
  Tdec_buf_out buffer_out; // output buffer
  buffer_out.reserve(2 * DEC_BUFFER_SIZE);
  std::size_t f{0};  // field
  std::size_t b{0};  // begin index of field in input buffer
  std::size_t i = b; // index in input buffer
  std::size_t constexpr N_FIELDS_SITE_DATA{9};

  // read data
  std::size_t bytes_read = fread(buffer_in.data(), 1, DEC_BUFFER_SIZE, stdin);
  std::vector<std::string> unique_fields{};
  unique_fields.reserve(16384);
  std::size_t remaining_bytes{0}; // bytes from previous buffer that weren't part of a complete field

  while (bytes_read != 0)
  {
    while (i < (bytes_read + remaining_bytes))
    {
      char const b_in = buffer_in[i];

      if (b_in != '\t' && b_in != '\n')
      {
        ++i;
        continue; // we are in a vcf field
      }

      if (f < N_FIELDS_SITE_DATA)
      {
        ++i; // adds '\t' or '\n'
        std::copy(buffer_in.begin() + b, buffer_in.begin() + i, std::back_inserter(buffer_out));
      }
      else
      {
        if (buffer_in[b] < ':')
        {
          unique_fields.emplace_back(&buffer_in[b], i - b); // add the unique field
          std::copy(buffer_in.begin() + b, buffer_in.begin() + (++i), std::back_inserter(buffer_out));
        }
        else
        {
          int32_t unique_index = ascii_cstring_to_int(&buffer_in[b], &buffer_in[i++]);
          assert(unique_index < static_cast<int32_t>(unique_fields.size()));
          std::string const & prior_field = unique_fields[unique_index];
          std::copy(prior_field.begin(), prior_field.end(), std::back_inserter(buffer_out));
          buffer_out.push_back(b_in);
        }
      }

      assert(b_in == buffer_in[i - 1]);
      b = i;

      if (b_in == '\n')
      {
        f = 0;
        unique_fields.clear();
      }
      else
      {
        ++f;
      }

      if (static_cast<long>(buffer_out.size()) >= DEC_BUFFER_SIZE)
      {
        fwrite(buffer_out.data(), 1, buffer_out.size(), stdout); // write buffer_out to stdout
        buffer_out.clear();
      }
    } // ends inner loop

    if (b == 0)
    {
      std::cerr << " ERROR: Encountered a field or line exceeding maximum buffer size of " << DEC_BUFFER_SIZE
                << std::endl;
      std::exit(1);
    }

    // write data to the beginning of the input buffer
    remaining_bytes = i - b;
    std::memmove(&buffer_in[0], &buffer_in[b], remaining_bytes);
    b = 0;
    i = remaining_bytes;
    bytes_read = fread(buffer_in.data() + remaining_bytes, 1, DEC_BUFFER_SIZE - remaining_bytes, stdin); // read more
  } /// ends outer loop

  fwrite(buffer_out.data(), 1, buffer_out.size(), stdout); // write to stdout
}

} // namespace popvcf
