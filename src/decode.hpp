#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <charconv>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include <parallel_hashmap/phmap.h>

#include "sequence_utils.hpp"

#include <popvcf/constants.hpp>

namespace popvcf
{
class DecodeData
{
public:
  std::size_t field{0};   //!< Current vcf field
  std::size_t in_size{0}; //!< Size of input buffer.
  std::size_t b{0};       //!< Field begin index in input buffer.
  std::size_t i{b};       //!< Curent index in input buffer
  bool header_line{true}; //!< True iff in header line
  bool in_region{true};   //!< True iff in region

  int64_t pos{-1};
  int64_t begin{-1};
  int64_t end{std::numeric_limits<int64_t>::max()};

  std::vector<uint32_t> prev_field2uid{};
  std::vector<std::string> prev_unique_fields{};
  phmap::flat_hash_map<std::string, uint32_t> prev_map_to_unique_fields{};

  std::vector<uint32_t> field2uid{};
  std::vector<std::string> unique_fields{};
  phmap::flat_hash_map<std::string, uint32_t> map_to_unique_fields{};

  inline void clear_line()
  {
    field = 0;
    std::swap(prev_field2uid, field2uid);
    std::swap(prev_unique_fields, unique_fields);
    std::swap(prev_map_to_unique_fields, map_to_unique_fields);
    field2uid.resize(0);
    unique_fields.resize(0);
    map_to_unique_fields.clear();
  }

  inline uint32_t num_unique_fields() const
  {
    return unique_fields.size();
  }
};

template <typename Tbuffer_in>
inline void set_input_size(Tbuffer_in & buffer_in, DecodeData & dd)
{
  dd.in_size = buffer_in.size();
}

template <>
inline void set_input_size(Tdec_array_buf & /*buffer_in*/, DecodeData & /*dd*/)
{
  // Do nothing.
  // NOTE: dd.in_size must be set prior to calling decode_buffer in arrays
}

//! Decodes an input buffer. Output is written in \a buffer_out .
template <bool is_region, typename Tbuffer_out, typename Tbuffer_in>
inline void decode_buffer(Tbuffer_out & buffer_out, Tbuffer_in & buffer_in, DecodeData & dd)
{
  set_input_size(buffer_in, dd);
  std::size_t constexpr N_FIELDS_SITE_DATA{9};

  // inner loop - Loops over each character in the input buffer
  while (dd.i < dd.in_size)
  {
    char const b_in = buffer_in[dd.i];

    if (b_in != '\t' && b_in != '\n')
    {
      ++dd.i; // we are in a vcf field
      continue;
    }

    if (dd.field == 0)
    {
      dd.header_line = buffer_in[dd.b] == '#'; // check if in header line
    }
    else if (dd.header_line == false && dd.field == 1)
    {
      std::from_chars(&buffer_in[dd.b], &buffer_in[dd.i], dd.pos); // get pos

      if constexpr (is_region)
      {
        dd.in_region = dd.pos >= dd.begin && dd.pos <= dd.end;
      }
    }

    if (dd.header_line || dd.field < N_FIELDS_SITE_DATA)
    {
      // write field without any encoding
      ++dd.i; // adds '\t' or '\n'

      if (!is_region || dd.in_region)
        buffer_out.insert(buffer_out.end(), &buffer_in[dd.b], &buffer_in[dd.i]);
    }
    else
    {
      long const field_idx = dd.field - N_FIELDS_SITE_DATA;
      assert(field_idx == static_cast<long>(dd.field2uid.size()));

      if (buffer_in[dd.b] == '$' || buffer_in[dd.b] == '&')
      {
        assert(field_idx < static_cast<long>(dd.prev_field2uid.size()));
        assert(dd.prev_field2uid[field_idx] < static_cast<long>(dd.prev_unique_fields.size()));

        std::string const & prior_field = dd.prev_unique_fields[dd.prev_field2uid[field_idx]];

        if (buffer_in[dd.b] == '$')
        {
          /* Unique field in this line. Same as field above. */
          dd.map_to_unique_fields.insert(std::pair<std::string, uint32_t>(prior_field, dd.num_unique_fields()));
          dd.field2uid.push_back(dd.num_unique_fields());
          dd.unique_fields.push_back(prior_field);
        }
        else
        {
          /* Duplicate field in this line. Same as field above. */
          auto find_it = dd.map_to_unique_fields.find(prior_field);
          assert(find_it != dd.map_to_unique_fields.end());
          dd.field2uid.push_back(find_it->second);
        }

        ++dd.i;

        if (!is_region || dd.in_region)
        {
          buffer_out.insert(buffer_out.end(), prior_field.begin(), prior_field.end());
          buffer_out.push_back(b_in);
        }
      }
      else if (buffer_in[dd.b] == '%')
      {
        // Unique field within the line but was seen in the previous line
        ++dd.b; // Get over '%'
        uint32_t const prev_unique_index = ascii_cstring_to_int(&buffer_in[dd.b], &buffer_in[dd.i++]);
        assert(prev_unique_index < dd.prev_unique_fields.size());
        std::string const & prior_field = dd.prev_unique_fields[prev_unique_index];

        dd.map_to_unique_fields.insert(std::pair<std::string, uint32_t>(prior_field, dd.num_unique_fields()));
        dd.field2uid.push_back(dd.num_unique_fields());
        dd.unique_fields.push_back(prior_field);

        if (!is_region || dd.in_region)
        {
          buffer_out.insert(buffer_out.end(), prior_field.begin(), prior_field.end());
          buffer_out.push_back(b_in);
        }
      }
      else if (buffer_in[dd.b] >= ':')
      {
        // same as earler field in the same line
        uint32_t const unique_index = ascii_cstring_to_int(&buffer_in[dd.b], &buffer_in[dd.i++]);
        assert(unique_index < dd.unique_fields.size());
        dd.field2uid.push_back(unique_index);
        std::string const & prior_field = dd.unique_fields[unique_index];

        if (!is_region || dd.in_region)
        {
          buffer_out.insert(buffer_out.end(), prior_field.begin(), prior_field.end());
          buffer_out.push_back(b_in);
        }
      }
      else
      {
        // add a new unique field and write field without any encoding
        auto insert_it = dd.map_to_unique_fields.insert(
          std::pair<std::string, uint32_t>(std::piecewise_construct,
                                           std::forward_as_tuple(&buffer_in[dd.b], dd.i - dd.b),
                                           std::forward_as_tuple(dd.num_unique_fields())));

        assert(insert_it.second == true);
        dd.field2uid.push_back(dd.num_unique_fields());
        dd.unique_fields.push_back(insert_it.first->first);
        ++dd.i;

        if (!is_region || dd.in_region)
          buffer_out.insert(buffer_out.end(), &buffer_in[dd.b], &buffer_in[dd.i]);
      }

      assert((field_idx + 1) == static_cast<long>(dd.field2uid.size()));
    }

    assert(b_in == buffer_in[dd.i - 1]);
    dd.b = dd.i;

    if (b_in == '\n')
      dd.clear_line();
    else
      ++dd.field;
  } // ends inner loop

  // write data to the beginning of the input buffer
  std::copy(&buffer_in[dd.b], &buffer_in[dd.i], &buffer_in[0]);
  dd.i = dd.i - dd.b;
  dd.b = 0;
  dd.in_size = dd.i;
  resize_input_buffer(buffer_in, dd.i);
}

//! Decode an encoded popVCF
void decode_file(std::string const & popvcf_fn, bool const is_bgzf_input);

//! Decode a region with a bgzf file and tabix index.
void decode_region(std::string const & popvcf_fn, std::string const & region);

} // namespace popvcf
