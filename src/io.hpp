#pragma once

#include <cstdio>
#include <iostream>
#include <string>

#include "htslib/bgzf.h"

class BGZF;

namespace popvcf
{
static inline BGZF * open_bgzf(std::string const & fn, std::string const & filemode)
{
  BGZF * in_bgzf = bgzf_open(fn.c_str(), filemode.c_str());

  if (in_bgzf == nullptr)
  {
    std::cerr << "[popvcf] ERROR: Opening bgzf file " << fn << std::endl;
    std::exit(1);
  }

  return in_bgzf;
}

static inline FILE * open_vcf(std::string const & fn, std::string const & filemode)
{
  if (fn == "-")
  {
    if (filemode == "r")
      return stdin;
    else
      return stdout;
  }

  FILE * in_vcf = fopen(fn.c_str(), filemode.c_str());

  if (in_vcf == nullptr)
  {
    std::cerr << "[popvcf] ERROR: Opening VCF file " << fn << std::endl;
    std::exit(1);
  }

  return in_vcf;
}

static inline void close_bgzf(BGZF * bgzf)
{
  int ret = bgzf_close(bgzf);

  if (ret != 0)
  {
    std::cerr << "[popvcf] ERROR: Failed closing bgzf file." << std::endl;
    std::exit(1);
  }
}
} // namespace popvcf
