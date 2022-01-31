#include "decode.hpp"

#include <algorithm> //std::copy
#include <array>     // std::array
#include <charconv>
#include <cstdio>   // std::stdin
#include <cstring>  // std::memmove
#include <iostream> // std::cerr
#include <string>   // std::string
#include <vector>   // std::vector

#include "io.hpp"
#include "sequence_utils.hpp" // ascii_cstring_to_int

#include "htslib/bgzf.h"

namespace popvcf
{
template <typename Tbuffer_out>
void decode_buffer(Tbuffer_out & buffer_out, Tdec_array_buf & buffer_in, DecodeData & dd)
{
  std::size_t constexpr N_FIELDS_SITE_DATA{9};

  // inner loop - Loops over each character in the input buffer
  while (dd.i < (dd.bytes_read + dd.remaining_bytes))
  {
    char const b_in = buffer_in[dd.i];

    if (b_in != '\t' && b_in != '\n')
    {
      // we are in a vcf field
      ++dd.i;
      continue;
    }

    if (dd.field == 0)
      dd.header_line = buffer_in[dd.b] == '#'; // check if in header line

    if (dd.header_line || dd.field < N_FIELDS_SITE_DATA)
    {
      // write field without any encoding
      ++dd.i; // adds '\t' or '\n'
      std::copy(buffer_in.begin() + dd.b, buffer_in.begin() + dd.i, std::back_inserter(buffer_out));
    }
    else
    {
      if (buffer_in[dd.b] == '$')
      {
	// write out a previous unique field in the same line
        uint32_t unique_index{0};
        std::from_chars(&buffer_in[dd.b + 1], &buffer_in[dd.i++], unique_index);
        std::string const & prior_field = dd.unique_fields[unique_index];
        std::copy(prior_field.begin(), prior_field.end(), std::back_inserter(buffer_out));
        buffer_out.push_back(b_in);
      }
      else if (buffer_in[dd.b] == '%')
      {
	// unique field within the line but was seen in the previous line
	uint32_t prev_unique_index{0};
        std::from_chars(&buffer_in[dd.b + 1], &buffer_in[dd.i++], prev_unique_index);
	assert(prev_unique_index < dd.prev_unique_fields.size());
        std::string const & prior_field = dd.prev_unique_fields[prev_unique_index];
        std::copy(prior_field.begin(), prior_field.end(), std::back_inserter(buffer_out));
	dd.unique_fields.emplace_back(prior_field);
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
    {
      dd.field = 0;
      std::swap(dd.prev_unique_fields, dd.unique_fields);
      dd.unique_fields.resize(0);
    }
    else
    {
      ++dd.field;
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
  std::copy(&buffer_in[dd.b], &buffer_in[dd.b + dd.remaining_bytes], &buffer_in[0]);
  dd.b = 0;
  dd.i = dd.remaining_bytes;
}

void decode_file(std::string const & input_fn, bool const is_bgzf_input)
{
  Tdec_array_buf buffer_in;     // input buffer
  std::vector<char> buffer_out; // output buffer
  DecodeData dd;                // data used to keep track of buffers while decoding

  /// Input streams
  BGZF * in_bgzf{nullptr};
  FILE * in_vcf{nullptr};

  /// Open input file based on options
  if (is_bgzf_input)
    in_bgzf = popvcf::open_bgzf(input_fn, "r");
  else
    in_vcf = popvcf::open_vcf(input_fn, "r");

  buffer_out.reserve(16 * DEC_BUFFER_SIZE);
  dd.unique_fields.reserve(32 * 1024);

  /// Read data
  if (is_bgzf_input)
    dd.bytes_read = bgzf_read(in_bgzf, buffer_in.data(), DEC_BUFFER_SIZE);
  else
    dd.bytes_read = fread(buffer_in.data(), 1, DEC_BUFFER_SIZE, in_vcf);

  /// Outer loop - loop while there is some data to decode from the input stream
  while (dd.bytes_read != 0)
  {
    decode_buffer(buffer_out, buffer_in, dd);

    /// Write buffer_out to stdout
    fwrite(buffer_out.data(), 1, buffer_out.size(), stdout);

    /// Clears output buffer, but does not deallocate
    buffer_out.resize(0);

    /// Read more data
    if (is_bgzf_input)
      dd.bytes_read = bgzf_read(in_bgzf, buffer_in.data() + dd.remaining_bytes, DEC_BUFFER_SIZE - dd.remaining_bytes);
    else
      dd.bytes_read = fread(buffer_in.data() + dd.remaining_bytes, 1, DEC_BUFFER_SIZE - dd.remaining_bytes, in_vcf);
  } /// ends outer loop

  assert(buffer_out.size() == 0);

  // fwrite(buffer_out.data(), 1, buffer_out.size(), stdout); // write to stdout

  if (is_bgzf_input)
    popvcf::close_bgzf(in_bgzf);
  else if (input_fn != "-")
    fclose(in_vcf);

  if (dd.remaining_bytes != 0)
  {
    std::cerr << "WARNING: After reading the input data we end inside of a field, the input may be truncated.\n"
              << "Remaining bytes=" << dd.remaining_bytes << " != 0" << std::endl;
  }
}

template void decode_buffer(std::vector<char> & buffer_out, Tdec_array_buf & buffer_in, DecodeData & dd);
template void decode_buffer(std::string & buffer_out, Tdec_array_buf & buffer_in, DecodeData & dd);

} // namespace popvcf
