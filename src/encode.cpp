#include "encode.hpp"

#include <array> // std::array
#include <charconv>
#include <iostream> // std::cerr
#include <string>   // std::string
#include <zlib.h>

#include <parallel_hashmap/phmap.h> // phmap::flat_hash_map

#include "io.hpp"
#include "sequence_utils.hpp" // int_to_ascii

#include "htslib/bgzf.h"

class BGZF;

namespace popvcf
{
void encode_file(std::string const & input_fn,
                 bool const is_bgzf_input,
                 std::string const & output_fn,
                 std::string const & output_mode,
                 bool const is_bgzf_output,
                 int const compression_threads)
{
  Tenc_array_buf buffer_in;     // input buffer
  std::vector<char> buffer_out; // output buffer
  EncodeData ed;                // encode data struct

  /// Open input file streams
  popvcf::bgzf_ptr in_bgzf(nullptr, popvcf::close_bgzf);   // bgzf input stream
  popvcf::file_ptr in_vcf(nullptr, popvcf::close_vcf_nop); // vcf input stream

  if (is_bgzf_input)
    in_bgzf = popvcf::open_bgzf(input_fn, "r");
  else
    in_vcf = popvcf::open_vcf(input_fn, "r");

  /// Open output file streams
  popvcf::bgzf_ptr out_bgzf(nullptr, popvcf::close_bgzf);   // bgzf output stream
  popvcf::file_ptr out_vcf(nullptr, popvcf::close_vcf_nop); // vcf output stream

  if (is_bgzf_output)
  {
    out_bgzf = popvcf::open_bgzf(output_fn.c_str(), output_mode.c_str());

    if (compression_threads > 1)
      bgzf_mt(out_bgzf.get(), compression_threads, 256);
  }
  else
  {
    out_vcf = popvcf::open_vcf(output_fn, output_mode);
  }

  /// Read first buffer of input data
  if (is_bgzf_input)
    ed.in_size = bgzf_read(in_bgzf.get(), buffer_in.data(), ENC_BUFFER_SIZE);
  else
    ed.in_size = fread(buffer_in.data(), 1, ENC_BUFFER_SIZE, in_vcf.get());

  long new_bytes = ed.in_size;

  // loop until all data has been read
  while (new_bytes != 0)
  {
    // encode the input buffer and write to output buffer
    encode_buffer(buffer_out, buffer_in, ed);

    // write output buffer
    if (out_bgzf != nullptr)
      popvcf::write_bgzf(out_bgzf.get(), buffer_out.data(), buffer_out.size());
    else
      fwrite(buffer_out.data(), 1, buffer_out.size(), out_vcf.get()); // write output buffer

    buffer_out.resize(0);
    new_bytes = -static_cast<long>(ed.in_size);

    // attempt to read more data from input
    if (is_bgzf_input)
      ed.in_size += bgzf_read(in_bgzf.get(), buffer_in.data() + ed.in_size, ENC_BUFFER_SIZE - ed.in_size);
    else
      ed.in_size += fread(buffer_in.data() + ed.in_size, 1, ENC_BUFFER_SIZE - ed.in_size, in_vcf.get());

    new_bytes += ed.in_size;
  }

  if (ed.in_size != 0)
  {
    std::cerr << "[popvcf] WARNING: Unexpected ending of the VCF data, possibly the file is truncated.\n";

    // write output buffer
    if (out_bgzf != nullptr)
      popvcf::write_bgzf(out_bgzf.get(), buffer_in.data(), ed.in_size);
    else
      fwrite(buffer_in.data(), 1, ed.in_size, out_vcf.get()); // write output buffer
  }
}

} // namespace popvcf
