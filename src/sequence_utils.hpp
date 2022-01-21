#pragma once

#include <cassert>
#include <string>
#include <string_view>
#include <vector>

namespace popvcf
{
int32_t constexpr CHAR_SET_SIZE = 69;
int32_t constexpr CHAR_SET_SIZE_2BYTES = CHAR_SET_SIZE * CHAR_SET_SIZE;
char constexpr CHAR_SET_MIN = ':';

inline char int_to_ascii(int32_t in)
{
  assert(in >= 0);
  assert(in < CHAR_SET_SIZE);

  return CHAR_SET_MIN + in;
}

inline int32_t ascii_to_int(char in)
{
  assert(in >= CHAR_SET_MIN);

  return static_cast<int32_t>(in) - CHAR_SET_MIN;
}

inline std::string int_to_ascii_string(int32_t in)
{
  assert(in >= 0);
  std::string str;

  while (in >= CHAR_SET_SIZE)
  {
    int32_t rem = in % CHAR_SET_SIZE;
    in = in / CHAR_SET_SIZE;
    str.push_back(int_to_ascii(rem));
  }

  assert(in < CHAR_SET_SIZE);
  str.push_back(int_to_ascii(in));
  return str;
}

inline int ascii_string_view_to_int(std::string_view in)
{
  int32_t const in_size = in.size();
  assert(in_size > 0);
  int32_t out{0};
  int32_t pow{1};

  for (int c{0}; c < in_size; ++c)
  {
    out += pow * ascii_to_int(in[c]);
    pow *= CHAR_SET_SIZE;
  }

  return out;
}

inline int32_t ascii_cstring_to_int(char const * b, char const * e)
{
  int32_t out{ascii_to_int(*b)};
  ++b;
  int32_t pow{CHAR_SET_SIZE};

  while (b != e)
  {
    out += pow * ascii_to_int(*b);
    pow *= CHAR_SET_SIZE;
    ++b;
  }

  return out;
}

template <typename Tstring>
std::vector<std::string_view> split_string(Tstring const & str, char const delimiter);

} // namespace popvcf
