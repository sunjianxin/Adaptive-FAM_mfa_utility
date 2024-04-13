project                     (example_MFA-DVR)
cmake_minimum_required      (VERSION 3.1)

option                      (wrapped_mpi            "MPI compiler wrapper requires no further MPI libs"     OFF     )
option                      (mfa_thread             "tbb, kokkos, sycl, or serial"                          serial  )
option                      (python                 "Build Python bindings"                                 OFF     )
option                      (MFA_INCLUDE_DIR        "path to MFA include directory"                                 )

set                         (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Default to Release
if                          (NOT CMAKE_BUILD_TYPE)
    set                     (CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
    set_property            (CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif                       (NOT CMAKE_BUILD_TYPE)

# OSX flags
if                          (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  add_definitions           (-DMAC_OSX)
  set                       (CMAKE_MACOSX_RPATH on)

    # --- following RPATH settings are for Sierra w/ Clang, hopefully they don't hurt other versions
    # ref: https://cmake.org/Wiki/CMake_RPATH_handling
    # use, i.e. don't skip, the full RPATH for the build tree
    set                     (CMAKE_SKIP_BUILD_RPATH            false)
    # when building, don't use the install RPATH already (but later on when installing)
    set                     (CMAKE_BUILD_WITH_INSTALL_RPATH    false)
    # set RPATH to install path
    set                     (CMAKE_INSTALL_RPATH               "${CMAKE_INSTALL_PREFIX}/lib")
    # add the automatically determined parts of the RPATH
    # which point to directories outside the build tree to the install RPATH
    set                     (CMAKE_INSTALL_RPATH_USE_LINK_PATH true)
    # the RPATH to be used when installing, but only if it's not a system directory
    list                    (FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
                                                               "${CMAKE_INSTALL_PREFIX}/lib"
                                                               isSystemDir)
    if                      ("${isSystemDir}" STREQUAL         "-1")
      set                   (CMAKE_INSTALL_RPATH               "${CMAKE_INSTALL_PREFIX}/lib")
    endif                   ()
endif                       (${CMAKE_SYSTEM_NAME} MATCHES      "Darwin")

message ("libraries:" ${libraries})

message ("---mpi include dir name:" ${MPI_INCLUDE_PATH})
message ("libraries:" ${libraries})
message ("mpi C libraries:" ${MPI_C_LIBRARIES})
message ("mpi CXX libraries:" ${MPI_CXX_LIBRARIES})

# MPI
 if                         (NOT wrapped_mpi)
     find_package           (MPI REQUIRED)
     set                    (libraries ${libraries} ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})
 endif                      ()

message ("---mpi include dir name:" ${MPI_INCLUDE_PATH})
message ("libraries:" ${libraries})
message ("mpi C libraries:" ${MPI_C_LIBRARIES})
message ("mpi CXX libraries:" ${MPI_CXX_LIBRARIES})

# MFA threading model
if                          (mfa_thread MATCHES "tbb")
    message                 ("Using TBB threading")
    add_definitions         (-DMFA_TBB)
    find_package            (TBB QUIET)
    find_path               (TBB_INCLUDE_DIR        tbb/tbb.h)
    find_library            (TBB_LIBRARY NAMES      tbb)
    if                      (TBB_INCLUDE_DIR AND TBB_LIBRARY)
        include_directories ("${TBB_INCLUDE_DIR}")
        set                 (libraries ${libraries} ${TBB_LIBRARY})
        message             ("Found TBB in ${TBB_INCLUDE_DIR} and in ${TBB_LIBRARY}")
    else                    (TBB_INCLUDE_DIR AND TBB_LIBRARY)
        message             (FATAL_ERROR "Could not find TBB")
    endif                   ()
else                        ()
    message                 ("Using no threading (serial)")
    add_definitions         (-DMFA_SERIAL)
endif                       ()

# Threads
find_package                (Threads)

message ("Threads libraries:" ${THREAD_INCLUDE_PATH})

# MFA
if                          (MFA_INCLUDE_DIR)
    include_directories     (SYSTEM ${MFA_INCLUDE_DIR})
    message                 ("MFA_INCLUDE_DIR = " ${MFA_INCLUDE_DIR})
    add_definitions         (-DMFA_NO_WEIGHTS)
    message                 ("Weights are disabled for MFA")
endif                       ()


# Include dirs
set                         (CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem")
include_directories         (${MFA_INCLUDE_DIR}/mfa
                             ${MFA_INCLUDE_DIR}/cppoptlib/include
                             ${MFA_INCLUDE_DIR}/diy/include
                             ${MFA_INCLUDE_DIR}/eigen3
                             SYSTEM
                             ${MPI_INCLUDE_PATH}
                             /usr/local/Cellar/boost/1.81.0_1/include
                             /Users/jianxinsun/research/intern/github/example_MFA-DVR/src
                            )

# Libraries
set                         (libraries
                            ${libraries}
                            ${CMAKE_DL_LIBS} 
                            ${CMAKE_THREAD_LIBS_INIT})

message ("-libraries:" ${libraries})
message ("-libraries:" ${CMAKE_DL_LIBS})
message ("-libraries:" ${CMAKE_THREAD_LIBS_INIT})
message ("CXX compiler:" ${CMAKE_CXX_COMPILER})
# Subdirectories
if                          (MFA_INCLUDE_DIR)
    add_subdirectory        (src)
else                        ()
    message                 (ERROR " Need MFA in order to build example_MFA-DVR")
endif                       ()

