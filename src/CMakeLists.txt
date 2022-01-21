# Update with "find src -name "*.?pp" | sort | awk '$1 !~ /main.cpp/{print "  "$1}'" in project root directory
set(popvcf_sources
  src/encode.cpp
  src/encode.hpp
  src/decode.cpp
  src/decode.hpp
  src/sequence_utils.cpp
  src/sequence_utils.hpp
  PARENT_SCOPE)
