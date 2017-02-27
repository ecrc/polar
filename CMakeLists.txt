cmake_minimum_required(VERSION 3.2.3)
# directly make an error if in-source build
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "In-source builds are not allowed.\n"
                        "Please create a build directory first and execute cmake configuration from "
                        "this directory. Example: mkdir build && cd build && cmake ..")
endif()

project(QDWH Fortran C)

#########################
### VARIABLES ###########
#########################
set( QDWH_LIBNAME qdwh )
option( QDWH_TESTING "Generates testing binaries" OFF)
option( BUILD_SHARED_LIBS "Generates shared libraries instead of static" OFF)
option( SHARED_NETLIBSCALAPACK "Build NETLIB-SCALAPACK as shared library." ${BUILD_SHARED_LIBS})
option( USE_PARALLEL_BLAS "If possible, prefer parallel blas libraries" ON )
set(EXTRA_LIBS "" CACHE STRING "Additional libraries that may be needed to compile/link" )

# Version
# Number 1.0.0 is temptative, may differ from actual release.
SET( QDWH_VERSION_MAJOR 1 )
SET( QDWH_VERSION_MINOR 0 )
SET( QDWH_VERSION_PATCH 0 )
set( QDWH_VERSION "${QDWH_VERSION_MAJOR}.${QDWH_VERSION_MINOR}.${QDWH_VERSION_PATCH}" )
#execute_process(COMMAND git rev-parse --short HEAD OUTPUT_VARIABLE QDWH_VERSION_RELEASE )


# Set the RPATH config
# --------------------
# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  FALSE)
# when building, use the install RPATH
# (automated test will need this)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
# the RPATH to be used when installing
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)


# Packaging (make package)
SET(CPACK_PACKAGE_VERSION ${QDWH_VERSION})
SET(CPACK_GENERATOR "TGZ")
INCLUDE(CPack)


# This project needs C99 standard to compile properly
# (could have used CMAKE_C_FLAGS instead)
#ADD_DEFINITIONS( -std=c99 )
#set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")

link_directories( ${CMAKE_INSTALL_PREFIX}/lib )
include_directories( ${CMAKE_SOURCE_DIR}/include )

## CMAKE MODULES :: ECRC
## REQUIRED FOR TESTS TO LINK LIBRARIES
if(NOT EXISTS "${PROJECT_SOURCE_DIR}/cmake_modules/ecrc/modules" )
    find_package(Git REQUIRED)
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule init WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} RESULT_VARIABLE _res_init OUTPUT_QUIET ERROR_QUIET)
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule update WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} RESULT_VARIABLE _res_update OUTPUT_QUIET ERROR_QUIET)
    if( ${_res_init} GREATER 0 OR ${_res_update} GREATER 0 )
        message(FATAL_ERROR "ECRC CMake modules were not found.\n"
                            "We tried: 'git submodule init && git submodule update' and resulted in error" )
    endif()
endif()
## ECRC INITIALIZATION
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake_modules/ecrc/modules")
set(ECRC_CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules/ecrc/modules )
include(EcrcInit)

# MPI
find_package( MPIEXT REQUIRED )
if( MPI_FOUND )
    # TODO check MPI_C_COMPILER
    include_directories( ${MPI_C_INCLUDE_PATH} )
    set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_LINK_FLAGS}" )
endif()

## BLAS
find_package(BLASEXT)
if(BLAS_FOUND)
# Use parallel blas
    if( USE_PARALLEL_BLAS AND NOT "${BLAS_PAR_LIBRARIES}" STREQUAL "" )
        message( STATUS "Using Parallel Blas" )
        set( BLAS_LIBRARIES "${BLAS_PAR_LIBRARIES}" )
    endif()
    if (BLAS_LIBRARY_DIRS)
        # the RPATH to be used when installing
        list( APPEND CMAKE_INSTALL_RPATH ${BLAS_LIBRARY_DIRS})
        link_directories( ${BLAS_LIBRARY_DIRS} )
    endif()
    if(BLAS_LINKER_FLAGS)
        list(APPEND CMAKE_EXE_LINKER_FLAGS ${BLAS_LINKER_FLAGS})
        list(APPEND CMAKE_SHARED_LINKER_FLAGS ${BLAS_LINKER_FLAGS})
        #set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}")
        list( REMOVE_DUPLICATES CMAKE_EXE_LINKER_FLAGS )
    endif()
else()
    message( FATAL_ERROR "BLAS library has not been found")
endif()

# LAPACK
#find_package(LAPACK REQUIRED)
find_package(LAPACKEXT REQUIRED)
if(LAPACK_FOUND)
    include_directories(${LAPACK_INCLUDE_DIRS})
    if(LAPACK_LIBRARY_DIRS)
        # the RPATH to be used when installing
        list(APPEND CMAKE_INSTALL_RPATH ${LAPACK_LIBRARY_DIRS})
        link_directories( ${LAPACK_LIBRARY_DIRS} )
    endif()
    if(LAPACK_LINKER_FLAGS)
        list( REMOVE_DUPLICATES LAPACK_LINKER_FLAGS )
        list( APPEND CMAKE_EXE_LINKER_FLAGS ${LAPACK_LINKER_FLAGS})
        list( REMOVE_DUPLICATES CMAKE_EXE_LINKER_FLAGS )
    endif()
endif()

# SCALAPACK
find_package( SCALAPACK )
if(SCALAPACK_FOUND)
    if (SCALAPACK_LIBRARY_DIRS)
        # the RPATH to be used when installing
        list(APPEND CMAKE_INSTALL_RPATH ${SCALAPACK_LIBRARY_DIRS})
        link_directories(${SCALAPACK_LIBRARY_DIRS})
    endif()
    if(SCALAPACK_LINKER_FLAGS)
        list( REMOVE_DUPLICATES SCALAPACK_LINKER_FLAGS )
        list( APPEND CMAKE_EXE_LINKER_FLAGS ${SCALAPACK_LINKER_FLAGS})
        list( REMOVE_DUPLICATES CMAKE_EXE_LINKER_FLAGS )
    endif()
endif()


include(ExternalProject)
# ELPA
find_package(ELPA)
if( ELPA_FOUND )
    include_directories( ${ELPA_INCLUDE_DIRS} )
    list(APPEND CMAKE_INSTALL_RPATH ${ELPA_LIBRARY_DIRS})
    link_directories( ${ELPA_LIBRARY_DIRS} )
else()
    set( ELPA_URL "http://elpa.mpcdf.mpg.de/html/Releases/2015.11.001/elpa-2015.11.001.tar.gz" )
    message( FATAL_ERROR "## ELPA not found. ##\n Download ELPA from ${ELPA_URL}, compile it and set PKG_CONFIG_PATH environment variable to the proper location" )
endif( ELPA_FOUND )


# SLTMG - look for sltmg (only used for tests )
if( QDWH_TESTING )
    set( SLTMG_LIBRARIES "" CACHE STRING "User provided sltmg library" )
    set( _notfound_msg_ "SLTMG not found, NETLIB-SCALAPACK will be downloaded in order to generate libsltmg.a\n If you want to use your own SLTMG, define SLTMG_LIBRARIES and make sure SCALAPACK is found." )
    if( SLTMG_LIBRARIES )
        set(CMAKE_REQUIRED_LIBRARIES "${SLTMG_LIBRARIES};${SCALAPACK_LIBRARIES}")
        check_fortran_function_exists("pdlatms" SLTMG_WORKS)
        set(CMAKE_REQUIRED_LIBRARIES)
        if( SLTMG_WORKS AND SCALAPACK_FOUND )
            message(STATUS "Found SLTMG : ${SLTMG_LIBRARIES}")
        else()
            message( STATUS "${_notfound_msg_}")
        endif()
    else()
        message( STATUS "${_notfound_msg_}")
    endif()
endif()

if( NOT SCALAPACK_FOUND OR ( QDWH_TESTING AND NOT SLTMG_WORKS ) )
    # Download NETLIB-SCALAPACK
    set( _build_netlib-scalapack )
    set( NETLIB_SCALAPACK_URL "http://www.netlib.org/scalapack/scalapack-2.0.2.tgz" CACHE STRING "URL to download NETLIB-SCALAPACK project" )
    if(NOT SCALAPACK_FOUND)
        set( SCALAPACK_LIBRARY_DIRS "${CMAKE_INSTALL_PREFIX}/lib")
        set( SCALAPACK_LIBRARIES "scalapack")
        set( _install_target "install" )
        set( _build_target "all" )
        set( _patch_file "scalapack-cmake-libsltmg.patch" )
        message( "## NETLIB-SCALAPACK will be downloaded from ${NETLIB_SCALAPACK_URL} and built along with libsltmg.a " )
    else()
        set( _install_target "install.sltmg" )
        set( _build_target "sltmg" )
        set( _patch_file "scalapack-cmake-libsltmg-only.patch" )
        message( "## NETLIB-SCALAPACK will be downloaded from ${NETLIB_SCALAPACK_URL} in order to generate libsltmg.a " )
    endif()

    # Use the built sltmg
    set( SLTMG_LIBRARIES "sltmg" )

    # build scalapack for sltmg and or elpa
    string(REPLACE ";" " " LAPACK_LIBRARIES_STR "${LAPACK_LIBRARIES}")
    ExternalProject_Add( NETLIB-SCALAPACK
                         PREFIX dependencies
                         URL ${NETLIB_SCALAPACK_URL}
                         #URL /tmp/scalapack-2.0.2.tgz
                         CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX} -DBUILD_SHARED_LIBS:BOOL=${SHARED_NETLIBSCALAPACK} -DLAPACK_LIBRARIES:STRING='${LAPACK_LIBRARIES_STR}'
                         PATCH_COMMAND patch -p1 -i ${CMAKE_SOURCE_DIR}/${_patch_file}
                         BUILD_COMMAND ${CMAKE_COMMAND} --build . --target ${_build_target} -- -j -l ${NUMBER_OF_CPU}
                         INSTALL_COMMAND ${CMAKE_COMMAND} --build . --target ${_install_target} --
                         #LOG_CONFIGURE 1
                         #LOG_BUILD 1
                         #LOG_INSTALL 1
    )

    if ( NOT ELPA_FOUND )
        # Build NETLIB-SCALAPACK before ELPA
        add_dependencies(ELPA NETLIB-SCALAPACK)
    endif()
endif()


# QDWH code
add_subdirectory(src)

# TESTS
if( QDWH_TESTING )

    ## LAPACKE
    find_package(LAPACKE COMPONENTS LAPACKEXT)
    if(LAPACK_FOUND AND LAPACK_LIBRARY_DIRS)
        # the RPATH to be used when installing
        list(APPEND CMAKE_INSTALL_RPATH ${LAPACK_LIBRARY_DIRS})
    else()
        message(FATAL_ERROR "A LAPACK library is required but has not been found")
    endif()
    if(LAPACKE_FOUND)
        if ( LAPACKE_INCLUDE_DIRS ) # Not always this variable is set, especially when is included in LAPACK
            include_directories(${LAPACKE_INCLUDE_DIRS})
        endif()
        if(LAPACKE_LIBRARY_DIRS)
            # the RPATH to be used when installing
            list(APPEND CMAKE_INSTALL_RPATH "${LAPACKE_LIBRARY_DIRS}")
        endif()
        if(LAPACKE_LINKER_FLAGS)
            list(APPEND CMAKE_EXE_LINKER_FLAGS ${LAPACKE_LINKER_FLAGS})
        endif()
    else()
        if(ECRC_VERBOSE_FIND_PACKAGE)
            if (LAPACKE_STANDALONE OR NOT LAPACKE_WORKS)
                if (NOT LAPACKE_lapacke.h_DIRS)
                    Print_Find_Header_Status(lapacke lapacke.h)
                endif ()
                if (NOT LAPACKE_lapacke_LIBRARY)
                    Print_Find_Library_Status(lapacke liblapacke)
                endif ()
            endif()
        else()
            message(WARNING "LAPACKE library has not been found and ECRC_VERBOSE_FIND_PACKAGE is set to OFF"
                            "Try to activate ECRC_VERBOSE_FIND_PACKAGE option (-DECRC_VERBOSE_FIND_PACKAGE=ON) to get some hints for the detection")
        endif()
        #message(FATAL_ERROR "A LAPACKE library is required but has not been found")
    endif()

    # BUILD TESTS
    enable_testing() # enables ctest
    add_definitions(-DHAVE_GETOPT_LONG)
    add_subdirectory(testing)
    add_subdirectory(timing)
endif()


# PKGCONFIG
set(prefix ${CMAKE_INSTALL_PREFIX})
set(libdir ${CMAKE_INSTALL_PREFIX}/lib)
set(PKG_CONFIG_DIR ${libdir}/pkgconfig)
configure_file(${CMAKE_SOURCE_DIR}/qdwh.pc.in ${CMAKE_CURRENT_BINARY_DIR}/qdwh.pc)
install( FILES ${CMAKE_BINARY_DIR}/qdwh.pc
         DESTINATION ${PKG_CONFIG_DIR}
       )
