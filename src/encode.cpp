#include "encode.hpp"

#include <array>    // std::array
#include <iostream> // std::cerr
#include <string>   // std::string
#include <zlib.h>

#include <parallel_hashmap/phmap.h> // phmap::flat_hash_map

#include "sequence_utils.hpp" // int_to_ascii

#include "htslib/bgzf.h"

namespace popvcf
{
using Tenc_map = phmap::flat_hash_map<std::string, int32_t>; //!< Hash table implementation to use when encoding

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

    if (ed.field < N_FIELDS_SITE_DATA)
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
        std::pair<std::string, int32_t>(std::piecewise_construct,
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
        int32_t char_val = insert_it.first->second; // character value of previous field, need to convert to string

        while (char_val >= CHAR_SET_SIZE)
        {
          auto rem = char_val % CHAR_SET_SIZE;
          char_val = char_val / CHAR_SET_SIZE;
          buffer_out[ed.o++] = int_to_ascii(rem);
        }

        assert(char_val < CHAR_SET_SIZE);
        buffer_out[ed.o++] = int_to_ascii(char_val);
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
                 std::string const & output_fn,
                 std::string const & output_mode,
                 bool const is_bgzf_output,
                 int const compression_threads)
{
  Tarray_buf buffer_in;  // input buffer
  Tarray_buf buffer_out; // output buffer
  EncodeData ed;         // encode data struct
  bool const is_stdin = input_fn == "-";

  // may be both gz file and normal plain text file
  gzFile fp = is_stdin ? nullptr : gzopen(input_fn.c_str(), "r");

  if (not stdin && fp == nullptr)
  {
    std::cerr << "[popvcf] ERROR: Could not open file " << input_fn << '\n';
    std::exit(1);
  }

  bool const is_stdout = output_fn == "-";
  BGZF * out_bgzf{nullptr};
  FILE * out{nullptr};

  // open output file based on options
  if (not is_bgzf_output)
  {
    if (not is_stdout)
    {
      out = fopen(output_fn.c_str(), output_mode.c_str());

      if (out == nullptr)
      {
        std::cerr << "[popvcf] ERROR: Opening file " << output_fn << std::endl;
        std::exit(1);
      }
    }
    else
    {
      out = stdout;
    }
  }
  else
  {
    out_bgzf = bgzf_open(output_fn.c_str(), output_mode.c_str());

    if (out_bgzf == nullptr)
    {
      std::cerr << "[popvcf] ERROR: Opening bgzf file " << output_fn << std::endl;
      std::exit(1);
    }

    if (compression_threads > 1)
      bgzf_mt(out_bgzf, compression_threads, 256);
  }

  // read first input data
  if (is_stdin)
    ed.bytes_read = fread(buffer_in.data(), 1, ENC_BUFFER_SIZE, stdin);
  else
    ed.bytes_read = gzread(fp, buffer_in.data(), ENC_BUFFER_SIZE);

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
      fwrite(buffer_out.data(), 1, ed.o, out); // write to stdout
    }

    // clear output index
    ed.o = 0;

    // attempt to read more data from input
    if (is_stdin)
    {
      ed.bytes_read = fread(buffer_in.data() + ed.remaining_bytes, 1, ENC_BUFFER_SIZE - ed.remaining_bytes, stdin);
    }
    else
    {
      ed.bytes_read = gzread(fp,                                    // file pointer
                             buffer_in.data() + ed.remaining_bytes, // buffer
                             ENC_BUFFER_SIZE - ed.remaining_bytes); // size
    }
  }

  if (not is_stdin)
    gzclose(fp);

  if (out_bgzf != nullptr)
  {
    int ret = bgzf_close(out_bgzf);

    if (ret != 0)
    {
      std::cerr << "[popvcf] ERROR: Failed closing bgzf file " << output_fn << std::endl;
      std::exit(1);
    }
  }
  else if (not stdout)
  {
    fclose(out);
  }
}

template void encode_buffer(std::string & buffer_out, Tarray_buf & buffer_in, EncodeData & ed);
template void encode_buffer(std::vector<char> & buffer_out, Tarray_buf & buffer_in, EncodeData & ed);
template void encode_buffer(Tarray_buf & buffer_out, Tarray_buf & buffer_in, EncodeData & ed);

} // namespace popvcf
