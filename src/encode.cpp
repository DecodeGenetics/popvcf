#include "encode.hpp"

#include <array>    // std::array
#include <iostream> // std::cerr
#include <string>   // std::string
#include <zlib.h>

#include <parallel_hashmap/phmap.h> // phmap::flat_hash_map

#include "sequence_utils.hpp" // int_to_ascii

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

void encode()
{
  Tarray_buf buffer_in;  // input buffer
  Tarray_buf buffer_out; // output buffer
  EncodeData ed;         // encode data struct

  // read the first buffer of input data from stdin
  ed.bytes_read = fread(buffer_in.data(), 1, ENC_BUFFER_SIZE, stdin);

  while (ed.bytes_read != 0)
  {
    // encode the input buffer and write to output buffer
    encode_buffer(buffer_out, buffer_in, ed);

    // write the output buffer to standard output
    fwrite(buffer_out.data(), 1, ed.o, stdout); // write to stdout
    ed.o = 0;

    // attempt to read more data from stdin
    ed.bytes_read = fread(buffer_in.data() + ed.remaining_bytes, // buffer
                          1,
                          ENC_BUFFER_SIZE - ed.remaining_bytes, // size
                          stdin);
  } /// ends outer loop
}

void encode_gzip_file(std::string const & gz_fn)
{
  Tarray_buf buffer_in;  // input buffer
  Tarray_buf buffer_out; // output buffer
  EncodeData ed;         // encode data struct

  gzFile fp = gzopen(gz_fn.c_str(), "r");

  if (fp == nullptr)
  {
    std::cerr << "ERROR: Could not open file " << gz_fn << '\n';
    std::exit(1);
  }

  ed.bytes_read = gzread(fp, buffer_in.data(), ENC_BUFFER_SIZE);

  while (ed.bytes_read != 0)
  {
    // encode the input buffer and write to output buffer
    encode_buffer(buffer_out, buffer_in, ed);

    // write the output buffer to standard output
    fwrite(buffer_out.data(), 1, ed.o, stdout); // write to stdout
    ed.o = 0;

    // attempt to read more data from stdin
    ed.bytes_read = gzread(fp,                                    // file pointer
                           buffer_in.data() + ed.remaining_bytes, // buffer
                           ENC_BUFFER_SIZE - ed.remaining_bytes); // size
  }
}

template void encode_buffer(std::string & buffer_out, Tarray_buf & buffer_in, EncodeData & ed);
template void encode_buffer(std::vector<char> & buffer_out, Tarray_buf & buffer_in, EncodeData & ed);
template void encode_buffer(Tarray_buf & buffer_out, Tarray_buf & buffer_in, EncodeData & ed);

} // namespace popvcf
