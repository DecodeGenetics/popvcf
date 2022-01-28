#pragma once
/*!
 * \file in.constants.hpp
 * \brief Global constants, macros and configurations set by CMake.
 */

#include <string.h> // strrchr

// clang-format off
// CMake variables
#define popvcf_VERSION_MAJOR @popvcf_VERSION_MAJOR@
#define popvcf_VERSION_MINOR @popvcf_VERSION_MINOR@
#define popvcf_VERSION_PATCH @popvcf_VERSION_PATCH@
#define popvcf_SOURCE_DIRECTORY "@PROJECT_SOURCE_DIR@"
#define popvcf_BINARY_DIRECTORY "@PROJECT_BINARY_DIR@"
#define GIT_BRANCH "@GIT_BRANCH@"
#define GIT_COMMIT_SHORT_HASH "@GIT_COMMIT_SHORT_HASH@"
#define GIT_COMMIT_LONG_HASH "@GIT_COMMIT_LONG_HASH@"
#define GIT_NUM_DIRTY_LINES "@GIT_NUM_DIRTY_LINES@"
// clang-format on

namespace popvcf
{
// Macros
#define S1_popvcf_internal__(x) #x
#define S2_popvcf_internal__(x) S1_popvcf_internal__(x)
#define _HERE_ (strrchr("/" __FILE__ ":" S2_popvcf_internal__(__LINE__), '/') + 1)

} // namespace popvcf
