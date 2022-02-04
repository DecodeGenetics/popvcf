#pragma once

#include <cstdio>
#include <iostream>
#include <memory>
#include <string>

#include "htslib/bgzf.h"
#include "htslib/hts.h"
#include "htslib/kseq.h"
#include "htslib/tbx.h"

class BGZF;

namespace popvcf
{
using file_ptr = std::unique_ptr<FILE, void (*)(FILE *)>;           //!< Type definition for a smart FILE pointer.
using bgzf_ptr = std::unique_ptr<BGZF, void (*)(BGZF *)>;           //!< Type definition for a smart BGZF pointer.
using hts_file_ptr = std::unique_ptr<htsFile, void (*)(htsFile *)>; //!< Type definition for a smart htsFile pointer.
using tbx_t_ptr = std::unique_ptr<tbx_t, void (*)(tbx_t *)>;        //!< Type definition for a smart tbx_t pointer.
using hts_itr_t_ptr = std::unique_ptr<hts_itr_t, void (*)(hts_itr_t *)>; //!< Type definition for a hts_itr_t pointer.

//! Closes a VCF file stream, i.e. stdout/stdin
inline void close_vcf_nop(FILE *)
{
}

//! Closes VCF file
inline void close_vcf(FILE * f)
{
  if (f != nullptr)
  {
    fclose(f);
  }
}

//! Opens a VCF file either from filename or stdout/stdin
inline file_ptr open_vcf(std::string const & fn, std::string const & filemode)
{
  if (fn == "-")
  {
    if (filemode == "r")
      return file_ptr(stdin, popvcf::close_vcf_nop);
    else
      return file_ptr(stdout, popvcf::close_vcf_nop);
  }

  file_ptr in_vcf(fopen(fn.c_str(), filemode.c_str()), popvcf::close_vcf);

  if (in_vcf == nullptr)
  {
    std::cerr << "[popvcf] ERROR: Opening VCF file " << fn << std::endl;
    std::exit(1);
  }

  return in_vcf;
}

inline void write_bgzf(BGZF * bgzf, const char * data, std::size_t const size)
{
  assert(bgzf != nullptr);
  std::size_t const written_bytes = bgzf_write(bgzf, data, size);

  if (written_bytes != size)
  {
    std::cerr << "[popvcf] WARNING: Problem writing bgzf data. " << written_bytes << " bytes written but expected "
              << size << " bytes.\n";
    std::exit(1);
  }
}

inline void close_bgzf(BGZF * bgzf)
{
  if (bgzf != nullptr)
  {
    if (bgzf_close(bgzf) != 0)
    {
      std::cerr << "[popvcf] ERROR: Failed closing bgzf file." << std::endl;
      std::exit(1);
    }
  }
}

inline bgzf_ptr open_bgzf(std::string const & fn, std::string const & filemode)
{
  bgzf_ptr in_bgzf(bgzf_open(fn.c_str(), filemode.c_str()), popvcf::close_bgzf);

  if (in_bgzf == nullptr)
  {
    std::cerr << "[popvcf] ERROR: Opening bgzf file " << fn << std::endl;
    std::exit(1);
  }

  return in_bgzf;
}

inline void close_hts_file(htsFile * f)
{
  if (f != nullptr)
  {
    if (hts_close(f) != 0)
    {
      std::cerr << "[popvcf] ERROR: Failed closing hts file." << std::endl;
      std::exit(1);
    }
  }
}

inline hts_file_ptr open_hts_file(const char * fn, const char * fm)
{
  hts_file_ptr ptr(hts_open(fn, fm), popvcf::close_hts_file);

  if (ptr == nullptr)
  {
    std::cerr << "ERROR: Could not open file " << fn << std::endl;
    std::exit(1);
  }

  return ptr;
}

inline void close_tbx_t(tbx_t * f)
{
  if (f != nullptr)
    tbx_destroy(f);
}

inline tbx_t_ptr open_tbx_t(const char * fn)
{
  tbx_t_ptr ptr(tbx_index_load(fn), popvcf::close_tbx_t);

  if (ptr == nullptr)
  {
    std::cerr << "[popvcf] ERROR: Could not open file " << fn << std::endl;
    std::exit(1);
  }

  return ptr;
}

inline void close_hts_itr_t(hts_itr_t * f)
{
  if (f != nullptr)
    tbx_itr_destroy(f);
}

inline hts_itr_t_ptr open_hts_itr_t(tbx_t * tbx, const char * region)
{
  hts_itr_t_ptr ptr(tbx_itr_querys(tbx, region), popvcf::close_hts_itr_t);

  if (ptr == nullptr)
    std::cerr << "[popvcf] WARNING: No records found in region " << region << "\n";

  return ptr;
}

inline void free_kstring_t(kstring_t * str)
{
  if (str->s != NULL)
    free(str->s);
}
} // namespace popvcf
