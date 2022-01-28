#include "decode.hpp"

#include <algorithm> //std::copy
#include <array>     // std::array
#include <cstring>   // std::memmove
#include <iostream>  // std::cerr
#include <string>    // std::string
#include <vector>    // std::vector
#include <zlib.h>

#include "sequence_utils.hpp" // int_to_ascii

namespace
{
} // namespace

namespace popvcf
{
// template <typename Tbuffer_out>
// void decode_buffer(Tbuffer_out & buffer_out, Tdec_array_buf & buffer_in, DecodeData & ed)
//{
//  // NOP
//}

void decode_file(std::string const & popvcf_fn)
{
  Tdec_array_buf buffer_in; // input buffer
  Tdec_buf_out buffer_out;  // output buffer
  DecodeData dd;
  bool const is_stdin = popvcf_fn == "-";
  gzFile in_fp = is_stdin ? nullptr : gzopen(popvcf_fn.c_str(), "r");

  if (not stdin && in_fp == nullptr)
  {
    std::cerr << "[popvcf] ERROR: Could not open file " << popvcf_fn << '\n';
    std::exit(1);
  }

  std::size_t constexpr N_FIELDS_SITE_DATA{9};
  buffer_out.reserve(2 * DEC_BUFFER_SIZE);

  // read data
  if (is_stdin)
    dd.bytes_read = fread(buffer_in.data(), 1, DEC_BUFFER_SIZE, stdin);
  else
    dd.bytes_read = gzread(in_fp, buffer_in.data(), DEC_BUFFER_SIZE);

  dd.unique_fields.reserve(16384);

  while (dd.bytes_read != 0)
  {
    while (dd.i < (dd.bytes_read + dd.remaining_bytes))
    {
      char const b_in = buffer_in[dd.i];

      if (b_in != '\t' && b_in != '\n')
      {
        ++dd.i;
        continue; // we are in a vcf field
      }

      if (dd.field < N_FIELDS_SITE_DATA)
      {
        ++dd.i; // adds '\t' or '\n'
        std::copy(buffer_in.begin() + dd.b, buffer_in.begin() + dd.i, std::back_inserter(buffer_out));
      }
      else
      {
        if (buffer_in[dd.b] < ':')
        {
          dd.unique_fields.emplace_back(&buffer_in[dd.b], dd.i - dd.b); // add the unique field
          std::copy(buffer_in.begin() + dd.b, buffer_in.begin() + (++dd.i), std::back_inserter(buffer_out));
        }
        else
        {
          int32_t unique_index = ascii_cstring_to_int(&buffer_in[dd.b], &buffer_in[dd.i++]);
          assert(unique_index < static_cast<int32_t>(dd.unique_fields.size()));
          std::string const & prior_field = dd.unique_fields[unique_index];
          std::copy(prior_field.begin(), prior_field.end(), std::back_inserter(buffer_out));
          buffer_out.push_back(b_in);
        }
      }

      assert(b_in == buffer_in[dd.i - 1]);
      dd.b = dd.i;

      if (b_in == '\n')
      {
        dd.field = 0;
        dd.unique_fields.clear();
      }
      else
      {
        ++dd.field;
      }

      if (static_cast<long>(buffer_out.size()) >= DEC_BUFFER_SIZE)
      {
        fwrite(buffer_out.data(), 1, buffer_out.size(), stdout); // write buffer_out to stdout
        buffer_out.clear();
      }
    } // ends inner loop

    if (dd.b == 0)
    {
      std::cerr << " ERROR: Encountered a field or line exceeding maximum buffer size of " << DEC_BUFFER_SIZE
                << std::endl;
      std::exit(1);
    }

    // write data to the beginning of the input buffer
    dd.remaining_bytes = dd.i - dd.b;
    std::memmove(&buffer_in[0], &buffer_in[dd.b], dd.remaining_bytes);
    dd.b = 0;
    dd.i = dd.remaining_bytes;

    if (is_stdin)
      dd.bytes_read = fread(buffer_in.data() + dd.remaining_bytes, 1, DEC_BUFFER_SIZE - dd.remaining_bytes, stdin);
    else
      dd.bytes_read = gzread(in_fp, buffer_in.data() + dd.remaining_bytes, DEC_BUFFER_SIZE - dd.remaining_bytes);
  } /// ends outer loop

  fwrite(buffer_out.data(), 1, buffer_out.size(), stdout); // write to stdout
}

} // namespace popvcf
