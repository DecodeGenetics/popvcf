#include "decode.hpp"

#include <algorithm> //std::copy
#include <array>     // std::array
#include <charconv>
#include <cstdio>   // std::stdin
#include <cstring>  // std::memmove
#include <iostream> // std::cerr
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

void decode_region(std::string const & popvcf_fn, std::string const & region)
{
  assert(region.size() > 0);
  std::string buffer_in; // input buffer
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
      auto ret = std::from_chars(&region[colon + 1], &region[region.size()], begin);

      if (ret.ec != std::errc())
        throw std::runtime_error("Could not parse region: " + region);

      end = begin;
    }
    else
    {
      auto ret_begin = std::from_chars(&region[colon + 1], &region[dash], begin);

      if (ret_begin.ec != std::errc())
        throw std::runtime_error("Could not parse region: " + region);

      auto ret_end = std::from_chars(&region[dash + 1], &region[region.size()], end);

      if (ret_end.ec != std::errc())
        throw std::runtime_error("Could not parse region: " + region);
    }

    dd.begin = begin;
    dd.end = end;
  }

  /// Input streams
  htsFile * in_bgzf = hts_open(popvcf_fn.c_str(), "r"); // open popvcf.gz

  if (in_bgzf == nullptr)
  {
    std::cerr << "ERROR: Could not open file " << popvcf_fn << std::endl;
    hts_close(in_bgzf);
    std::exit(1);
  }

  tbx_t * in_tbx = tbx_index_load(popvcf_fn.c_str()); // open popvcf.gz.tbi

  if (in_tbx == nullptr)
  {
    std::cerr << "ERROR: Could not open file " << popvcf_fn << std::endl;
    std::exit(1);
  }

  /// Copy the header lines
  kstring_t str = {0, 0, 0};

  while (hts_getline(in_bgzf, KS_SEP_LINE, &str) >= 0)
  {
    if (!str.l || str.s[0] != in_tbx->conf.meta_char)
    {
      break;
    }

    std::copy(str.s, str.s + str.l, std::back_inserter(buffer_out));
    buffer_out.push_back('\n');
  }

  /// Write buffer_out to stdout
  fwrite(buffer_out.data(), 1, buffer_out.size(), stdout);

  /// Clears output buffer, but does not deallocate
  buffer_out.resize(0);

  /// Query the region
  std::string safe_region = chrom;

  if (begin >= 0)
  {
    safe_region.push_back(':');
    safe_region.append(std::to_string((begin / 10000) * 10000));
    safe_region.push_back('-');
    safe_region.append(std::to_string(end));
  }

  hts_itr_t * in_it = tbx_itr_querys(in_tbx, safe_region.c_str());

  if (in_it == nullptr)
  {
    std::cerr << "ERROR: Could query region " << safe_region << std::endl;
    hts_close(in_bgzf);
    tbx_destroy(in_tbx);
    std::exit(1);
  }

  int ret = tbx_itr_next(in_bgzf, in_tbx, in_it, &str);

  while (ret > 0)
  {
    buffer_in.resize(dd.remaining_bytes);
    std::copy(str.s, str.s + str.l, std::back_inserter(buffer_in));
    buffer_in.push_back('\n');
    dd.bytes_read = str.l + 1;

    decode_buffer(buffer_out, buffer_in, dd);

    /// Write buffer_out to stdout
    fwrite(buffer_out.data(), 1, buffer_out.size(), stdout);

    /// Clears output buffer, but does not deallocate
    buffer_out.resize(0);

    /// Check if end position has been passed
    if (dd.pos > dd.end)
      break;

    /// Read more data
    ret = tbx_itr_next(in_bgzf, in_tbx, in_it, &str);
  }

  hts_close(in_bgzf);
  tbx_destroy(in_tbx);
  tbx_itr_destroy(in_it);
}

template void decode_buffer(std::vector<char> & buffer_out, Tdec_array_buf & buffer_in, DecodeData & dd);
template void decode_buffer(std::string & buffer_out, Tdec_array_buf & buffer_in, DecodeData & dd);

} // namespace popvcf
