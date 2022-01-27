#include <cstdio>
#include <cstdlib>
#include <iostream>

#include <paw/parser.hpp>

#include "decode.hpp"
#include "encode.hpp"

#include <popvcf/constants.hpp>

namespace popvcf
{
int subcmd_encode(paw::Parser & parser)
{
  std::string vcf_fn{};

  try
  {
    parser.parse_positional_argument(vcf_fn,
                                     "VCF",
                                     "Encode this VCF (or VCF.gz). If not set, read VCF from standard input.");

    parser.finalize();
  }
  catch (paw::exception::missing_positional_argument &)
  {
  }

  if (vcf_fn.empty() || vcf_fn == "-")
    encode();
  else
    encode_gzip_file(vcf_fn);

  return 0;
}

int subcmd_decode(paw::Parser & /*parser*/)
{
  // std::string popvcf_fn{};
  // parser.parse_positional_argument(popvcf_fn, "popVCF", "Decode this popVCF. Use '-' for standard input.");

  decode();
  return 0;
}

} // namespace popvcf

int main(int argc, char ** argv)
{
#ifndef NDEBUG
  std::ios_base::sync_with_stdio(false);
#endif // NDEBUG
  paw::Parser parser(argc, argv);
  parser.set_name("popVCF");
  parser.set_version(popvcf_VERSION_MAJOR, popvcf_VERSION_MINOR, popvcf_VERSION_PATCH);
  int ret{0};

  try
  {
    std::string subcmd{};

    parser.add_subcommand("encode", "Encode a VCF into a popVCF.");
    parser.add_subcommand("decode", "Decode a popVCF into a VCF.");

    parser.parse_subcommand(subcmd);

    if (subcmd == "encode")
    {
      ret = popvcf::subcmd_encode(parser);
    }
    else if (subcmd == "decode")
    {
      ret = popvcf::subcmd_decode(parser);
    }
    else if (subcmd.size() == 0)
    {
      parser.finalize();
      ret = 0;
    }
    else
    {
      parser.finalize();
      ret = 1;
    }
  }
  catch (paw::exception::help const & e)
  {
    std::cout << e.what();
    return 0;
  }
  catch (std::exception const & e)
  {
    std::cerr << e.what();
    return 1;
  }

  return ret;
}
