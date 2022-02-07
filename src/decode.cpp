#include "decode.hpp"

#include <algorithm> //std::copy
#include <array>     // std::array
#include <charconv>
#include <cstdio>   // std::stdin
#include <cstring>  // std::memmove
#include <iostream> // std::cerr
#include <memory>
#include <stdexcept>
#include <string> // std::string
#include <vector> // std::vector

#include "io.hpp"
#include "sequence_utils.hpp" // ascii_cstring_to_int

#include "htslib/bgzf.h"
#include "htslib/hts.h"
#include "htslib/kseq.h"
#include "htslib/tbx.h"

namespace popvcf
{
void decode_file(std::string const & input_fn, bool const is_bgzf_input)
{
  Tdec_array_buf buffer_in;     // input buffer
  std::vector<char> buffer_out; // output buffer
  DecodeData dd;                // data used to keep track of buffers while decoding

  /// Input streams
  popvcf::bgzf_ptr in_bgzf(nullptr, popvcf::close_bgzf);
  popvcf::file_ptr in_vcf(nullptr, popvcf::close_vcf_nop);

  /// Open input file based on options
  if (is_bgzf_input)
    in_bgzf = popvcf::open_bgzf(input_fn, "r");
  else
    in_vcf = popvcf::open_vcf(input_fn, "r");

  buffer_out.reserve(16 * DEC_BUFFER_SIZE);

  /// Read first batch of data
  if (is_bgzf_input)
    dd.in_size = bgzf_read(in_bgzf.get(), buffer_in.data(), DEC_BUFFER_SIZE);
  else
    dd.in_size = fread(buffer_in.data(), 1, DEC_BUFFER_SIZE, in_vcf.get());

  long new_bytes = dd.in_size;

  /// Outer loop - loop while there is some data to decode from the input stream
  while (new_bytes != 0)
  {
    decode_buffer</*in_region=*/false>(buffer_out, buffer_in, dd);

    /// Write buffer_out to stdout
    fwrite(buffer_out.data(), 1, buffer_out.size(), stdout);
    buffer_out.resize(0); // Clears output buffer, but does not deallocate
    new_bytes = -static_cast<long>(dd.in_size);

    /// Read more data
    if (is_bgzf_input)
      dd.in_size += bgzf_read(in_bgzf.get(), buffer_in.data() + dd.in_size, DEC_BUFFER_SIZE - dd.in_size);
    else
      dd.in_size += fread(buffer_in.data() + dd.in_size, 1, DEC_BUFFER_SIZE - dd.in_size, in_vcf.get());

    new_bytes += dd.in_size;
  } /// ends outer loop

  assert(buffer_out.size() == 0);

  if (dd.in_size != 0)
  {
    std::cerr << "[popvcf] WARNING: Unexpected ending of the VCF data, possibly the file is truncated.\n";

    // write output buffer
    fwrite(buffer_in.data(), 1, dd.in_size, stdout); // write output buffer
  }
}

void decode_region(std::string const & popvcf_fn, std::string const & region)
{
  assert(region.size() > 0);
  std::vector<char> buffer_in; // input buffer
  buffer_in.reserve(DEC_BUFFER_SIZE);
  std::vector<char> buffer_out; // output buffer
  DecodeData dd;                // data used to keep track of buffers while decoding

  /// parse region
  std::string chrom;
  long begin{-1};
  long end{std::numeric_limits<long>::max()};

  if (auto colon = region.find(':'); colon == std::string::npos)
  {
    chrom = region;
  }
  else
  {
    chrom = region.substr(0, colon);

    if (auto dash = region.find('-', colon + 1); dash == std::string::npos)
    {
      auto ret = std::from_chars(region.data() + colon + 1, region.data() + region.size(), begin);

      if (ret.ec != std::errc())
        throw std::runtime_error("Could not parse region: " + region);

      end = begin;
    }
    else
    {
      auto ret_begin = std::from_chars(region.data() + colon + 1, region.data() + dash, begin);

      if (ret_begin.ec != std::errc())
        throw std::runtime_error("Could not parse region: " + region);

      auto ret_end = std::from_chars(region.data() + dash + 1, region.data() + region.size(), end);

      if (ret_end.ec != std::errc())
        throw std::runtime_error("Could not parse region: " + region);
    }

    dd.begin = begin;
    dd.end = end;
  }

  /// Determine the region to query
  std::string safe_region = chrom;

  if (begin >= 0)
  {
    safe_region.push_back(':');
    safe_region.append(std::to_string(std::max(1l, (begin / 10000l) * 10000l)));
    safe_region.push_back('-');
    safe_region.append(std::to_string(end));
  }

  /// Input streams
  popvcf::hts_file_ptr in_bgzf = popvcf::open_hts_file(popvcf_fn.c_str(), "r");            // open popvcf.gz
  popvcf::tbx_t_ptr in_tbx = popvcf::open_tbx_t(popvcf_fn.c_str());                        // open popvcf.gz.tbi
  popvcf::hts_itr_t_ptr in_it = popvcf::open_hts_itr_t(in_tbx.get(), safe_region.c_str()); // query region

  /// Write the header lines
  kstring_t str = {0, 0, 0};

  while (hts_getline(in_bgzf.get(), KS_SEP_LINE, &str) >= 0)
  {
    if (!str.l || str.s[0] != in_tbx->conf.meta_char)
      break;

    fwrite(str.s, 1, str.l, stdout);
    fputs("\n", stdout);
  }

  // return here, after writing header, if there are no records in the region
  if (in_it == nullptr)
  {
    free(str.s);
    return;
  }

  int ret = tbx_itr_next(in_bgzf.get(), in_tbx.get(), in_it.get(), &str);

  while (ret > 0)
  {
    buffer_in.insert(buffer_in.end(), str.s, str.s + str.l);
    buffer_in.push_back('\n');

    decode_buffer</*in_region=*/true>(buffer_out, buffer_in, dd);

    /// Write buffer_out to stdout
    fwrite(buffer_out.data(), 1, buffer_out.size(), stdout);

    /// Clears output buffer, but does not deallocate
    buffer_out.resize(0);

    /// Check if end position has been passed
    if (dd.pos > dd.end)
      break;

    /// Read more data
    ret = tbx_itr_next(in_bgzf.get(), in_tbx.get(), in_it.get(), &str);
  }

  free(str.s);
}

} // namespace popvcf
