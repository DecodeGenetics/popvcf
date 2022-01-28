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
  std::string vcf_fn{"-"};
  std::string input_type{"v"};
  std::string output_fn{"-"};
  std::string output_mode{"w"};
  std::string output_type{"v"};
  int output_compression_level{-1};
  int compression_threads{1};

  try
  {
    parser.parse_positional_argument(vcf_fn,
                                     "VCF",
                                     "Encode this VCF (or VCF.gz). If not set, read VCF from standard input.");

    parser.parse_option(input_type, 'I', "input-type", "v|z", "Input type. v uncompressed VCF, z bgzipped VCF.");
    parser.parse_option(output_fn, 'o', "output", "VCF.gz", "Output filename.");
    parser.parse_option(output_compression_level,
                        'l',
                        "output-compress-level",
                        "INT",
                        "Output file compression level.");
    parser.parse_option(compression_threads, '@', "threads", "INT", "Output file compression level.");
    parser.parse_option(output_type, 'O', "output-type", "v|z", "Output type. v uncompressed VCF, z bgzipped VCF.");

    parser.finalize();
  }
  catch (paw::exception::missing_positional_argument &)
  {
    output_fn = "-";
  }

  if (output_compression_level >= 0)
    output_mode += std::to_string(std::min(9, output_compression_level));

  encode_file(vcf_fn, input_type == "z", output_fn, output_mode, output_type == "z", compression_threads);
  return 0;
}

int subcmd_decode(paw::Parser & parser)
{
  std::string popvcf_fn{};
  std::string input_type{"v"};

  try
  {
    parser.parse_option(input_type, 'I', "input-type", "v|z", "Input type. v uncompressed VCF, z bgzipped VCF.");
    parser.parse_positional_argument(popvcf_fn, "popVCF", "Decode this popVCF. Use '-' for standard input.");
    parser.finalize();
  }
  catch (paw::exception::missing_positional_argument &)
  {
    popvcf_fn = "-";
  }

  long const n = popvcf_fn.size();

  if (n > 3 && popvcf_fn[n - 2] == 'g' && popvcf_fn[n - 1] == 'z')
    input_type = "z";

  decode_file(popvcf_fn, input_type == "z");
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
