
cmake_minimum_required(VERSION 3.5.1)
project(BParser)

message(STATUS "CMakeLists.txt - BParser")

# CLANG
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC") #clang-cl
  message(STATUS "Compiler: clang-cl")
  string(CONCAT CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} " /std:c++17 /Wall "
  "/WX " #missing "/clang:-Wno-*" or "-Xclang -Wno-*"
  "-DBOOST_PHOENIX_NO_VARIADIC_EXPRESSION") 
  #https://clang.llvm.org/docs/UsersManual.html#the-clang-option
  #https://learn.microsoft.com/en-us/cpp/build/reference/compiler-options?view=msvc-170#find-a-compiler-option

  string(CONCAT CMAKE_CXX_FLAGS_DEBUG  ${CMAKE_CXX_FLAGS_DEBUG} " -gcodeview /Zi ")# /MT /fsanitize=address /clang:-shared-libasan") #linker issues with asan => unsupported

else()
  message(STATUS "Compiler: clang/gnu") #other compilers will likely fail
  #set(CMAKE_CXX_FLAGS     "-std=c++14 -finline-hint-functions -pedantic-errors -Werror=pedantic -Wall -Wextra -Werror -Wno-long-long -Wno-strict-aliasing -DBOOST_PHOENIX_NO_VARIADIC_EXPRESSION")
  #set(CMAKE_CXX_FLAGS     "-std=c++14 -fPIC -pedantic-errors -Werror=pedantic -Wall -Wextra -Werror -Wno-long-long -Wno-strict-aliasing -DBOOST_PHOENIX_NO_VARIADIC_EXPRESSION")
  string(CONCAT CMAKE_CXX_FLAGS  ${CMAKE_CXX_FLAGS} " -std=c++17 -finline-functions -pedantic-errors -Werror=pedantic -Wall -Wextra "
  "-Werror -Wno-unused-function -Wno-unused-variable -Wno-long-long -Wno-strict-aliasing "
  "-Wno-microsoft-cpp-macro -Wno-microsoft-enum-value "
  "-DBOOST_PHOENIX_NO_VARIADIC_EXPRESSION ")

  if(ON) # use sanitizer
      # valgrind setting
      string(CONCAT CMAKE_CXX_FLAGS_DEBUG  ${CMAKE_CXX_FLAGS_DEBUG}  " -g -Og -DDEBUG  -fsanitize=address -fno-omit-frame-pointer")  #-static-libasan
      if(WIN32) 
        string(CONCAT CMAKE_CXX_FLAGS_DEBUG  ${CMAKE_CXX_FLAGS_DEBUG} " -shared-libasan")
      endif()
  endif()

  if(FLOW_BUILD_TYPE STREQUAL "debug")
    set(CMAKE_CXX_FLAGS_DEBUG   "-g -DDEBUG -fno-omit-frame-pointer")  #
  endif()
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3")


# TODO: Try -Winline to be detect possibly non inlined functions.


#ASAN_OPT         := -g -fsanitize=address -static-libasan -fno-omit-frame-pointer
#MSAN_OPT         := -g -fsanitize=memory    -fno-omit-frame-pointer
#LSAN_OPT         := -g -fsanitize=leak      -fno-omit-frame-pointer
#USAN_OPT         := -g -fsanitize=undefined -fno-omit-frame-pointer


# Set a default build type if none was specified
set(default_build_type "Release")
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
  set(default_build_type "Debug")
endif()
 
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo" "Sanitizer")
endif()

# TODO: modify makefile to call cmake
# - cmake_debug
# - cmake_release
# - all
# - clean

######################################################################### find libraries
#
# find boost or use our own
#
# - Try to find Boost include dirs and libraries
#
# NOTE: Take note of the Boost_ADDITIONAL_VERSIONS variable below.
# Due to Boost naming conventions and limitations in CMake this find
# module is NOT future safe with respect to Boost version numbers,
# and may break.
#
# Input variables:
#
#   Boost_USE_MULTITHREADED      Can be set to OFF to use the non-multithreaded
#   Boost_USE_STATIC_LIBS        Can be set to ON to force the use of the static
#   Boost_NO_SYSTEM_PATHS        Set to TRUE to suppress searching in system (or other locations outside of BOOST_ROOT or BOOST_INCLUDEDIR).
#   Boost_USE_STATIC_RUNTIME     If enabled, searches for boost libraries linked against a static C++ standard library
#   Boost_USE_DEBUG_PYTHON       
#   Boost_USE_STLPORT            If enabled, searches for boost libraries compiled against the STLPort standard
#   Boost_USE_STLPORT_DEPRECATED_NATIVE_IOSTREAMS    If enabled, searches for boost libraries compiled against the deprecated STLPort "native iostreams" feature ('n' ABI tag).
#   Boost_ADDITIONAL_VERSIONS    A list of version numbers to use for searching
#   Boost_DEBUG                  Set this to TRUE to enable debugging output
#   Boost_DETAILED_FAILURE_MSG   FindBoost doesn't output detailed information
#   Boost_COMPILER               Set this to the compiler suffix used by Boost
#   Boost_THREADAPI                When building boost.thread, sometimes the name of the library contains an additional "pthread" or "win32"
#   Boost_REALPATH               Resolves symbolic links for discovered boost libraries
#
#
#   BOOST_ROOT or BOOSTROOT      The preferred installation prefix for searching for
#                                back on the system paths, set Boost_NO_SYSTEM_PATHS
#   BOOST_INCLUDEDIR             Set this to the include directory of Boost, if the
#   BOOST_LIBRARYDIR             Set this to the lib directory of Boost, if the
#
# Variables defined by this module:
#
#   Boost_FOUND                         System has Boost, this means the include dir was
#   Boost_INCLUDE_DIRS                  Boost include directories: not cached
#   Boost_INCLUDE_DIR                   This is almost the same as above, but this one is cached and may be modified by advanced users
#   Boost_LIBRARIES                     Link to these to use the Boost libraries that you
#   Boost_LIBRARY_DIRS                  The path to where the Boost library files are.
#   Boost_VERSION                       The version number of the boost libraries that
#   Boost_LIB_VERSION                   The version number in filename form as
#   Boost_MAJOR_VERSION                 major version number of boost
#   Boost_MINOR_VERSION                 minor version number of boost
#   Boost_SUBMINOR_VERSION              subminor version number of boost
#   Boost_LIB_DIAGNOSTIC_DEFINITIONS    [WIN32 Only] You can call
#   Boost_${COMPONENT}_FOUND            True IF the Boost library "component" was found.
#   Boost_${COMPONENT}_LIBRARY          Contains the libraries for the specified Boost

message(STATUS "=======================================================")
message(STATUS "====== BOOST ==========================================")
message(STATUS "=======================================================")
set(Boost_USE_STATIC_LIBS       ON)
set(Boost_USE_MULTITHREADED     OFF)
set(Boost_USE_STATIC_RUNTIME    OFF)
set(EXTERNAL_BOOST_DIR  "${EXTERNAL_PROJECT_DIR}/boost_build")  
 
# set(Boost_DEBUG TRUE)
# For some mysterious reason, the second find_package fails to find
# components of boost in our own build, if we hint BOOST_ROOT of our build to the
# first find_package.
#
# The root of the error is in FindBoost, when it calls find_library for individual components.
# Success scenario:
# 1) First find_package didn't found Boost (in system direcotries)
# 2) we install our Boost
# 3) the second find_package calls find_library to find boost components and succed
#
# Fail scenario:
# 1) First find_package didn't found Boost (we set BOOST_ROOT to the empty target directory of the external project)
# 2) we install our Boost
# 3) the second find_package calls find_library to find boost components and fail
#    the parameters of the find_library are exactly same as in the previous case.
#

# If no explicit setting from config.cmake ....
if (NOT BOOST_ROOT)
    set(BOOST_ROOT $ENV{BOOST_ROOT})
endif()    

# If no system setting ...
if (NOT BOOST_ROOT)
    set(BOOST_ROOT ${EXTERNAL_BOOST_DIR})
endif()

if (NOT Boost_FOUND) 
  # we try to erease from the cache everything that can cause failure of the second find_package
  # however just unsetting Boost_INCLUDE_DIR should be enough
 
  unset(Boost_FOUND CACHE) 
  unset(Boost_INCLUDE_DIR CACHE)        # overwrite cache    
  unset(Boost_LIBRARY_DIRS CACHE)
  unset(Boost_LIB_VERSION CACHE)
  unset(Boost_VERSION CACHE)
  
  set(BOOST_ROOT ${EXTERNAL_BOOST_DIR})    
  if (NOT EXTERNAL_PROJECT_DIR)
   unset(BOOST_ROOT)
  endif()
  find_package( Boost 1.58.0 REQUIRED)
endif()

message(STATUS "-------------------------------------------------------")
message(STATUS "Boost_FOUND = ${Boost_FOUND}")
message(STATUS "BOOST_ROOT = ${BOOST_ROOT}")
message(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}")
message(STATUS "Boost_LIBRARY_DIRS = ${Boost_LIBRARY_DIRS}")
message(STATUS "Boost_INCLUDE_DIR = ${Boost_INCLUDE_DIR}")
message(STATUS "=======================================================\n\n")

message(STATUS "VCL2_INCLUDE_DIR = ${CMAKE_CURRENT_SOURCE_DIR}/third_party/VCL_v2")

set(BPARSER_INCLUDES   ${CMAKE_CURRENT_SOURCE_DIR}/include     ${Boost_INCLUDE_DIR}     ${CMAKE_CURRENT_SOURCE_DIR}/third_party/VCL_v2)
if(NOT PROJECT_IS_TOP_LEVEL)
    set(BPARSER_INCLUDES ${BPARSER_INCLUDES} PARENT_SCOPE)
endif()


include_directories(
    ${BPARSER_INCLUDES}
)

# add tager for libbparser

set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/include/processor_SSE.cc -ffast-math -mfpmath=sse -msse4)
set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/include/processor_AVX2.cc -ffast-math -mfpmath=sse -mfma -mavx2)
set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/include/processor_AVX512.cc -ffast-math -mfpmath=sse -mfma -mavx512f)

add_library(bparser SHARED
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/VCL_v2/instrset_detect.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/grammar.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/include/processor_SSE.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/include/processor_AVX2.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/include/processor_AVX512.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/include/processor_double.cc
)







######################################
# Test targets


#################################
# MACRO define_test
#
# Creates build targets for particular source file containing one test.
#
# For given parameter 'class_name', it assumes class_name_test.cpp source
# and creates:
# target 'class_name_test' which builds corresponding binary and run the test
# target 'class_name_valgrind' which builds the binary and runs it under valgrind
# use libraries in ${super_test_libs}
#
# We want to run tests from source dir, but since there is no portable way
# how to make link from source dir to Makefile generated in build dir
# I let main makefile to make such links.
macro(define_test src_name)
  
  set(test_source "${CMAKE_CURRENT_SOURCE_DIR}/test/${src_name}.cc")
  set(test_binary "${src_name}_bin")
  set(test_name "${src_name}")
  set(test_valgrind "v${src_name}")
  #message(STATUS "adding test: ${test_binary} cxx flags: ${CMAKE_CXX_FLAGS}")

  # executable rule
  add_executable(${test_binary} EXCLUDE_FROM_ALL  ${test_source} )
  if(CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC") #CL interface (Clang-cl)
    set_target_properties(${test_binary} PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} /WX-")
  else()
    set_target_properties(${test_binary} PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -Wno-error")
  endif()
  set(test_dependencies ${ARGN})
  if (test_dependencies) 
        add_dependencies(${test_binary} ${test_dependencies}) # all additional macro arguments are the test dependencies
        target_link_libraries(${test_binary} ${test_dependencies})
  endif()      
  
  # rule to run the test
  add_custom_target(${test_name}
      COMMAND "$<TARGET_FILE:${test_binary}>"
      DEPENDS ${test_binary}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  # rule to run the test under valgrind
  add_custom_target(${test_valgrind}
      COMMAND valgrind --num-callers=50 "${CMAKE_CURRENT_BINARY_DIR}/${test_binary}"
      DEPENDS ${test_binary}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endmacro()


define_test(test_parser bparser)
define_test(test_array)
define_test(test_grammar bparser)
define_test(test_processor)
define_test(test_speed bparser)  
define_test(test_simd)
