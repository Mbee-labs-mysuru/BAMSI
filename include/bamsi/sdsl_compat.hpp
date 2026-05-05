#pragma once

// Only bother on C++17+ where invoke_result exists.
#include <type_traits>

#if __cplusplus >= 201703L

// On libc++ (AppleClang) also try to re-enable removed C++17 features.
#  if defined(__clang__) && defined(_LIBCPP_VERSION)
#    ifndef _LIBCPP_ENABLE_CXX17_REMOVED_FEATURES
#      define _LIBCPP_ENABLE_CXX17_REMOVED_FEATURES
#    endif
#    ifndef _LIBCPP_ENABLE_CXX17_REMOVED_RESULT_OF
#      define _LIBCPP_ENABLE_CXX17_REMOVED_RESULT_OF
#    endif
#  endif

// If the library advertises is_invocable / invoke_result, redirect result_of.
#  if defined(__cpp_lib_is_invocable) && __cpp_lib_is_invocable >= 201703
     // This macro turns `std::result_of<...>` into `std::invoke_result<...>`.
#    ifndef SDSL_BAMSI_RESULT_OF_SHIM
#      define SDSL_BAMSI_RESULT_OF_SHIM
#      define result_of invoke_result
#    endif
#  endif

#endif
