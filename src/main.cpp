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
  std::string input_type{"g"};
  std::string output_fn{"-"};
  std::string output_mode{"w"};
  std::string output_type{"v"};
  int output_compress_level{-1};
  int compression_threads{1};

  try
  {
    parser.parse_positional_argument(vcf_fn,
                                     "VCF",
                                     "Encode this VCF (or VCF.gz). If not set, read VCF from standard input.");

    parser.parse_option(compression_threads,
                        '@',
                        "threads",
                        "Number of output file compression threads (only used if output type is \"z\").",
                        "NUM");

    parser.parse_option(input_type,
                        'I',
                        "input-type",
                        "Input type. v uncompressed VCF, z bgzipped VCF, g guess based on filename.",
                        "v|z|g");

    parser.parse_option(output_fn,
                        'o',
                        "output",
                        "Output will be written to this path. If '-', then write instead to standard output.",
                        "output.vcf[.gz]");

    parser.parse_option(output_compress_level, 'l', "output-compress-level", "Output file compression level.", "LEVEL");

    parser.parse_option(output_type, 'O', "output-type", "Output type. v uncompressed VCF, z bgzipped VCF.", "v|z");
    parser.finalize();
  }
  catch (paw::exception::missing_positional_argument &)
  {
    output_fn = "-";
  }

  if (output_compress_level >= 0)
    output_mode += std::to_string(std::min(9, output_compress_level));

  long const n = vcf_fn.size();

  if (n > 3 && vcf_fn[n - 2] == 'g' && vcf_fn[n - 1] == 'z')
    input_type = "z";

  encode_file(vcf_fn, input_type == "z", output_fn, output_mode, output_type == "z", compression_threads);
  return 0;
}

int subcmd_decode(paw::Parser & parser)
{
  std::string popvcf_fn{};
  std::string input_type{"g"};
  std::string region{};

  try
  {
    parser.parse_option(input_type,
                        'I',
                        "input-type",
                        "Input type. v uncompressed VCF, z bgzipped VCF, g guess based on filename.",
                        "v|z|g");
    parser.parse_option(region, 'r', "region", "Fetch region/interval to decode. Requires .tbi index.", "chrN:A-B");
    parser.parse_positional_argument(popvcf_fn, "popVCF", "Decode this popVCF. Use '-' for standard input.");
    parser.finalize();
  }
  catch (paw::exception::missing_positional_argument &)
  {
    popvcf_fn = "-";
  }

  long const n = popvcf_fn.size();

  if (input_type == "g" && n > 3 && popvcf_fn[n - 2] == 'g' && popvcf_fn[n - 1] == 'z')
    input_type = "z";

  if (region.empty())
    decode_file(popvcf_fn, input_type == "z");
  else
    decode_region(popvcf_fn, region);

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
