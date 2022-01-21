#include <cstdio>
#include <cstdlib>
#include <iostream>

#include <paw/parser.hpp>

#include "decode.hpp"
#include "encode.hpp"

#include <popvcf/constants.hpp>

namespace popvcf
{
int subcmd_encode(paw::Parser & /*parser*/)
{
  // std::string vcf_fn{};
  // parser.parse_positional_argument(vcf_fn, "VCF", "Encode this VCF. Use '-' for standard input.");

  encode();
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
    parser.add_subcommand("getline", "getline on input (for benchmarking)");
    parser.add_subcommand("fread", "fread on input (for benchmarking)");

    parser.parse_subcommand(subcmd);

    if (subcmd == "encode")
    {
      ret = popvcf::subcmd_encode(parser);
    }
    else if (subcmd == "decode")
    {
      ret = popvcf::subcmd_decode(parser);
    }
    else if (subcmd == "getline")
    {
      for (std::string line; std::getline(std::cin, line);)
        std::cout << line << '\n';

      ret = 0;
    }
    else if (subcmd == "fread")
    {
      long constexpr size{65536};
      char * buffer;
      buffer = static_cast<char *>(malloc(size));

      while (fread(buffer, 1, size, stdin) != 0)
      {
        std::cout << buffer;
      }

      free(buffer);
      ret = 0;
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
