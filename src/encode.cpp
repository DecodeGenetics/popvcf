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
                 int const compression_threads,
                 bool const no_previous_line)
{
  Tarray_buf buffer_in;  // input buffer
  Tarray_buf buffer_out; // output buffer
  EncodeData ed;         // encode data struct
  ed.no_previous_line = no_previous_line;

  // Open input file streams
  BGZF * in_bgzf{nullptr}; // bgzf input stream
  FILE * in_vcf{nullptr};  // vcf input stream

  if (is_bgzf_input)
    in_bgzf = popvcf::open_bgzf(input_fn, "r");
  else
    in_vcf = popvcf::open_vcf(input_fn, "r");

  // Open output file streams
  BGZF * out_bgzf{nullptr}; // bgzf output stream
  FILE * out_vcf{nullptr};  // vcf output stream

  if (is_bgzf_output)
  {
    out_bgzf = popvcf::open_bgzf(output_fn.c_str(), output_mode.c_str());

    if (compression_threads > 1)
      bgzf_mt(out_bgzf, compression_threads, 256);
  }
  else
  {
    out_vcf = popvcf::open_vcf(output_fn, output_mode);
  }

  /// Read first buffer of input data
  if (is_bgzf_input)
    ed.bytes_read = bgzf_read(in_bgzf, buffer_in.data(), ENC_BUFFER_SIZE);
  else
    ed.bytes_read = fread(buffer_in.data(), 1, ENC_BUFFER_SIZE, in_vcf);

  // loop until all data has been read
  while (ed.bytes_read != 0)
  {
    // encode the input buffer and write to output buffer
    encode_buffer(buffer_out, buffer_in, ed);

    // write output buffer
    if (out_bgzf != nullptr)
    {
      long const written_bytes = bgzf_write(out_bgzf, buffer_out.data(), ed.o);

      if (written_bytes != static_cast<long>(ed.o))
      {
        std::cerr << "[popvcf] WARNING: Problem writing bgzf data to " << output_fn << " " << written_bytes
                  << " bytes written but expected " << ed.o << " bytes." << std::endl;
      }
    }
    else
    {
      fwrite(buffer_out.data(), 1, ed.o, out_vcf); // write output buffer
    }

    ed.o = 0; // clear output index

    // attempt to read more data from input
    if (is_bgzf_input)
      ed.bytes_read = bgzf_read(in_bgzf, buffer_in.data() + ed.remaining_bytes, ENC_BUFFER_SIZE - ed.remaining_bytes);
    else
      ed.bytes_read = fread(buffer_in.data() + ed.remaining_bytes, 1, ENC_BUFFER_SIZE - ed.remaining_bytes, in_vcf);
  }

  /// Close input streams
  if (in_bgzf != nullptr)
    popvcf::close_bgzf(in_bgzf);
  else if (output_fn != "-")
    fclose(in_vcf);

  /// Close output streams
  if (out_bgzf != nullptr)
    popvcf::close_bgzf(out_bgzf);
  else if (output_fn != "-")
    fclose(out_vcf);

  if (ed.remaining_bytes != 0)
  {
    std::cerr << "WARNING: After reading the input data we end inside of a field, the input may be truncated.\n"
              << "Remaining bytes=" << ed.remaining_bytes << " != 0" << std::endl;
  }
}

/*
template void encode_buffer(std::string & buffer_out, Tarray_buf & buffer_in, EncodeData & ed);
template void encode_buffer(std::vector<char> & buffer_out, Tarray_buf & buffer_in, EncodeData & ed);
template void encode_buffer(Tarray_buf & buffer_out, Tarray_buf & buffer_in, EncodeData & ed);
template void encode_buffer(std::string & buffer_out, std::string & buffer_in, EncodeData & ed);
template void encode_buffer(std::vector<char> & buffer_out, std::string & buffer_in, EncodeData & ed);
template void encode_buffer(std::string & buffer_out, std::vector<char> & buffer_in, EncodeData & ed);
template void encode_buffer(std::vector<char> & buffer_out, std::vector<char> & buffer_in, EncodeData & ed);
*/

} // namespace popvcf
