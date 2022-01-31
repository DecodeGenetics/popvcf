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
template <typename Tbuffer_out>
void reserve_space(Tbuffer_out & buffer)
{
  buffer.resize(ENC_BUFFER_SIZE);
}

//! Specialization for array buffers
template <>
void reserve_space(Tarray_buf & /*buffer*/)
{
  // NOP
}

template <typename Tbuffer_out>
void encode_buffer(Tbuffer_out & buffer_out, Tarray_buf & buffer_in, EncodeData & ed)
{
  reserve_space(buffer_out);
  std::size_t constexpr N_FIELDS_SITE_DATA{9}; // how many fields of the VCF contains site data

  while (ed.i < (ed.bytes_read + ed.remaining_bytes))
  {
    char const b_in = buffer_in[ed.i];

    if (b_in != '\t' && b_in != '\n')
    {
      ++ed.i;
      continue; // we are in a vcf field
    }

    if (ed.field == 0)
      ed.header_line = buffer_in[ed.b] == '#'; // check if in header line

    if (ed.header_line || ed.field < N_FIELDS_SITE_DATA /*|| (ed.i - ed.b) <= 5*/)
    {
      ++ed.i; // adds '\t' or '\n'
      std::copy(&buffer_in[ed.b], &buffer_in[ed.i], &buffer_out[ed.o]);
      ed.o += (ed.i - ed.b);
    }
    else
    {
      assert(buffer_in[ed.b] >= '!');
      assert(buffer_in[ed.b] <= '9');

      auto insert_it = ed.map_to_unique_fields.insert(
        std::pair<std::string, uint32_t>(std::piecewise_construct,
                                         std::forward_as_tuple(&buffer_in[ed.b], ed.i - ed.b),
                                         std::forward_as_tuple(ed.num_unique_fields)));

      if (insert_it.second == true)
      {
        ++ed.num_unique_fields; // unique field
        ++ed.i;                 // adds '\t' or '\n'
        std::copy(&buffer_in[ed.b], &buffer_in[ed.i], &buffer_out[ed.o]);
        ed.o += (ed.i - ed.b);
      }
      else
      {
        // field identical to prior field
        buffer_out[ed.o++] = '$';
        auto ret = std::to_chars(&buffer_out[ed.o], &buffer_out[ed.o + 9], insert_it.first->second);

#ifndef NDEBUG
        if (ret.ec == std::errc::value_too_large)
        {
          std::cerr << "ERROR: Index value too large\n";
          std::exit(1);
        }
#endif // NDEBUG
        ed.o += (ret.ptr - &buffer_out[ed.o]);
        buffer_out[ed.o++] = buffer_in[ed.i++];
      }
    }

    assert(b_in == buffer_in[ed.i - 1]); // i should have been already incremented here
    ed.b = ed.i;                         // set begin index of next field

    // check if we need to clear line or increment field
    if (b_in == '\n')
      ed.clear_line();
    else
      ++ed.field;
  } // ends inner loop

  if (ed.b == 0)
  {
    std::cerr << " ERROR: Encountered a field or line exceeding maximum buffer size of " << ENC_BUFFER_SIZE
              << std::endl;
    std::exit(1);
  }

  // copy the remaining data to the beginning of the input buffer
  ed.remaining_bytes = ed.i - ed.b;
  std::copy(&buffer_in[ed.b], &buffer_in[ed.b + ed.remaining_bytes], &buffer_in[0]);
  ed.b = 0;
  ed.i = ed.remaining_bytes;
}

void encode_file(std::string const & input_fn,
                 bool const is_bgzf_input,
                 std::string const & output_fn,
                 std::string const & output_mode,
                 bool const is_bgzf_output,
                 int const compression_threads)
{
  Tarray_buf buffer_in;  // input buffer
  Tarray_buf buffer_out; // output buffer
  EncodeData ed;         // encode data struct

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

template void encode_buffer(std::string & buffer_out, Tarray_buf & buffer_in, EncodeData & ed);
template void encode_buffer(std::vector<char> & buffer_out, Tarray_buf & buffer_in, EncodeData & ed);
template void encode_buffer(Tarray_buf & buffer_out, Tarray_buf & buffer_in, EncodeData & ed);

} // namespace popvcf
