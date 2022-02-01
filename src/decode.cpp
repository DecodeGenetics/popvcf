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
