######################################
# ---------------------------------- #
# -------- COMPILER VERSION -------- #
# ---------------------------------- #
######################################

cmake_minimum_required(VERSION 3.11)
project(receptor_detector VERSION 0.0.0)

message(STATUS "CMAKE_CXX_COMPILER_ID: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION: ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_SYSTEM: ${CMAKE_SYSTEM}")

# Compiler options.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (dev)
    message(STATUS "---- Building with debug options ----")

    # Memory debug
    set(CMAKE_CXX_FLAGS "-fexceptions -fsanitize=address")
    set(ASAN_OPTIONS=check_initialization_order=1)
    set(ASAN_OPTIONS=detect_leaks=1)

    add_definitions(-ggdb3 -O0 -Wall)       # Debugging + No optimization

else()
    message(STATUS "---- Building with optimization ----")

    # Standard compilation
    add_definitions(-O3 -Wall)              # Much optimization
endif()



#########################################
# ------------------------------------- #
# -------- SOURCES AND HEADERS -------- #
# ------------------------------------- #
#########################################

# Include header files
include_directories(
        "inc"
)

# Define our shared library sources. NOT test/executables.
set(SOURCES
        src/BinaryIO.cpp
        src/BinarySequence.cpp
        src/Hasher2.cpp
        src/IncrementalIdMap.cpp
        src/misc.cpp
        src/MurmurHash2.cpp
        src/MurmurHash3.cpp
        src/Sequence.cpp
        src/SvgPlot.cpp
        src/Timer.cpp
        src/TsvReader.cpp
        )


project(receptor_detector)
add_library(receptor_detector STATIC ${SOURCES})

# To make sure the library is named Name.so,
# get rid of the "lib" prefix.
set_target_properties(receptor_detector PROPERTIES PREFIX "")

# Eliminate an extraneous -D during compilation.
set_target_properties(receptor_detector PROPERTIES DEFINE_SYMBOL "")


##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ----------------- edlib ------------------ #
# ------------------------------------------ #
##############################################

# Need to explicitly enable ExternalProject functionality
include(ExternalProject)

# Download or update library as an external project
ExternalProject_Add(project_edlib
        GIT_REPOSITORY https://github.com/Martinsos/edlib.git
        PREFIX ${CMAKE_SOURCE_DIR}/external/edlib
        #        DOWNLOAD_COMMAND ""
        #        UPDATE_COMMAND ""
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external/edlib/ -DCMAKE_INSTALL_LIBDIR=${CMAKE_SOURCE_DIR}/external/edlib/lib
        STAMP_DIR ${CMAKE_SOURCE_DIR}/external/edlib/stamp
        #        SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/edlib/src
        BINARY_DIR ${CMAKE_SOURCE_DIR}/external/edlib/build
        INSTALL_DIR ${CMAKE_SOURCE_DIR}/external/edlib/
        INSTALL_COMMAND make install
        )

# Define INSTALL_DIR as the install directory for external library
ExternalProject_Get_Property(project_edlib INSTALL_DIR)

# Create new library for external project (so it can be linked with main library)
add_library(edlib STATIC IMPORTED)
set_property(TARGET edlib
        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/lib/libedlib.a)

# Define main library as dependent on the downloaded project (transitively)
add_dependencies(edlib
        project_edlib)

add_dependencies(receptor_detector edlib)

# Ensure that main library has access to primary dependencies' and secondary dependencies' headers
include_directories(external/edlib/include/)

# Add runtime path for main library so it can access external library
#set_property(TARGET GetBlunted PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/bin/ogdf")

message(STATUS "INSTALL_DIR: ${INSTALL_DIR}")


##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ---------------- sparsepp ---------------- #
# ------------------------------------------ #
##############################################

include(ExternalProject)

ExternalProject_Add(project_sparsepp
        GIT_REPOSITORY https://github.com/greg7mdp/sparsepp.git
        TIMEOUT 5
        PREFIX "${CMAKE_CURRENT_BINARY_DIR}"
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND "" # Disable install step, is a header only lib!
        )

add_dependencies(receptor_detector project_sparsepp)

# Specify include dir
ExternalProject_Get_Property(project_sparsepp SOURCE_DIR)
include_directories(${SOURCE_DIR}/)
include_directories(${SOURCE_DIR}/sparsepp)

############################################
# ---------------------------------------- #
# -------- Generating executables -------- #
# ---------------------------------------- #
############################################

#
##set(Boost_DEBUG 1)
#find_package( Boost 1.60 COMPONENTS program_options REQUIRED )
#include_directories( ${Boost_INCLUDE_DIR} )
#
#

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# -------- TESTS --------

set(TESTS
        test_binary_sequence
        test_binary_sequence_performance
        test_edlib
        test_kmer_classifier
        )

foreach(FILENAME_PREFIX ${TESTS})
    add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
    target_link_libraries(${FILENAME_PREFIX}
            receptor_detector
            Threads::Threads
            edlib
            )

endforeach()


# -------- EXECUTABLES --------

set(EXECUTABLES
        compute_minhash
        compute_rle_minhash
        )

foreach(FILENAME_PREFIX ${EXECUTABLES})
    add_executable(${FILENAME_PREFIX} src/executable/${FILENAME_PREFIX}.cpp)
    target_link_libraries(${FILENAME_PREFIX}
            receptor_detector
            Threads::Threads
            edlib
            )

endforeach()

#set_target_properties(receptor_detector PROPERTIES LINK_FLAGS "-static" )
#SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++ -static-libgcc")


# -------- final steps --------

# Where to install
set(BINARY_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)

## The shared library goes to the bin directory.
#install(TARGETS
#        detector
#        receptor_detector
#        DESTINATION ${BINARY_INSTALL_DIR})

# Some potentially helpful messages
message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")
