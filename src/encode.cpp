#include "encode.hpp"

#include <array>    // std::array
#include <iostream> // std::cerr
#include <string>   // std::string

#include <parallel_hashmap/phmap.h> // phmap::flat_hash_map

#include "sequence_utils.hpp" // int_to_ascii

namespace popvcf
{
using Tenc_map = phmap::flat_hash_map<std::string, int32_t>; //!< Hash table implementation to use when encoding
using Tbuf = std::array<char, ENC_BUFFER_SIZE>;              //!< Buffer type

void encode()
{
  Tbuf buffer_in;    // input buffer
  Tbuf buffer_out;   // output buffer
  std::size_t o{0};  // output index
  std::size_t f{0};  // field
  std::size_t b{0};  // begin index in buffer_in
  std::size_t i = b; // index in buffer_in
  std::size_t constexpr N_FIELDS_SITE_DATA{9};

  // read data
  std::size_t bytes_read = fread(buffer_in.data(), 1, ENC_BUFFER_SIZE, stdin);
  int32_t num_unique_fields{0};
  Tenc_map map_to_unique_fields{};
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
        std::memcpy(&buffer_out[o], &buffer_in[b], i - b);
        o += (i - b);
      }
      else
      {
        auto insert_it =
          map_to_unique_fields.insert(std::pair<std::string, int32_t>(std::piecewise_construct,
                                                                      std::forward_as_tuple(&buffer_in[b], i - b),
                                                                      std::forward_as_tuple(num_unique_fields)));

        if (insert_it.second == true)
        {
          ++num_unique_fields; // unique field
          ++i;                 // adds '\t' or '\n'
          std::memcpy(&buffer_out[o], &buffer_in[b], i - b);
          o += (i - b);
        }
        else
        {
          // field identical to prior field
          int32_t char_val = insert_it.first->second; // character value of previous field, need to convert to string

          while (char_val >= CHAR_SET_SIZE)
          {
            auto rem = char_val % CHAR_SET_SIZE;
            char_val = char_val / CHAR_SET_SIZE;
            buffer_out[o++] = int_to_ascii(rem);
          }

          assert(char_val < CHAR_SET_SIZE);
          buffer_out[o++] = int_to_ascii(char_val);
          buffer_out[o++] = buffer_in[i++];
        }
      }

      assert(b_in == buffer_in[i - 1]); // i should have been incremented here
      b = i;

      if (b_in == '\n')
      {
        f = 0;
        map_to_unique_fields.clear();
        num_unique_fields = 0;
      }
      else
      {
        ++f;
      }
    } // ends inner loop

    if (b == 0)
    {
      std::cerr << " ERROR: Encountered a field or line exceeding maximum buffer size of " << ENC_BUFFER_SIZE
                << std::endl;
      std::exit(1);
    }

    fwrite(buffer_out.data(), 1, o, stdout); // write to stdout
    o = 0;

    // write data to the beginning of the input buffer
    remaining_bytes = i - b;
    std::memmove(&buffer_in[0], &buffer_in[b], remaining_bytes);
    b = 0;
    i = remaining_bytes;
    bytes_read = fread(buffer_in.data() + remaining_bytes, 1, ENC_BUFFER_SIZE - remaining_bytes, stdin); // read more
  } /// ends outer loop
}

} // namespace popvcf
