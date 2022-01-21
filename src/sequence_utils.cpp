#include "sequence_utils.hpp"

#include <algorithm>   // std::find
#include <cstdint>     // int32_t
#include <string>      // std::string
#include <string_view> // std::string_view
#include <vector>      // std::vector

namespace popvcf
{
template <typename Tstring>
std::vector<std::string_view> split_string(Tstring const & str, char const delimiter)
{
  std::vector<std::string_view> output;
  std::string_view strv(str);
  auto first = strv.cbegin();

  while (first != strv.cend())
  {
    auto const second = std::find(first, strv.cend(), delimiter);

    if (first != second)
    {
      std::size_t const pos = std::distance(strv.cbegin(), first);
      output.emplace_back(strv.substr(pos, second - first));
    }

    if (second == strv.cend())
      break;

    first = std::next(second);
  }

  return output;
}

template std::vector<std::string_view> split_string(std::string const & str, char const delimiter);
template std::vector<std::string_view> split_string(std::string_view const & str, char const delimiter);

} // namespace popvcf
