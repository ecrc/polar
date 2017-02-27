# Defines the source files of the tests. Each file generates 1 test
file(GLOB tests_files "*.c")

include_directories( ${CMAKE_SOURCE_DIR}/include )
link_directories( ${CMAKE_BINARY_DIR}/lib )

# Uses RUNPATH instead of RPATH
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--enable-new-dtags")
list( APPEND CMAKE_EXE_LINKER_FLAGS "-Wl,--enable-new-dtags")
# Transform the list into a space separate string
string(REPLACE ";" " " CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

foreach( test_src ${tests_files} )
    get_filename_component( test_exe ${test_src}  NAME_WE )
    if( "${test_exe}" STREQUAL "timing" )
        # skip timing.c
        continue()
    endif()

    add_executable( ${test_exe} timing.c ${test_src} )
    add_dependencies( ${test_exe} ${QDWH_LIBNAME} )
    if( _build_netlib-scalapack )
        add_dependencies( ${test_exe} NETLIB-SCALAPACK )
    endif()

    target_link_libraries( ${test_exe} ${QDWH_LIBNAME}
                                       ${MPI_C_LIBRARIES}
                                       ${ELPA_LIBRARIES}
                                       ${SLTMG_LIBRARIES}
                                       ${SCALAPACK_LIBRARIES}
                                       ${LAPACKE_LIBRARIES}
                                       ${LAPACK_LIBRARIES}
                                       ${BLAS_LIBRARIES}
                                       ${EXTRA_LIBS} )
    install( TARGETS ${test_exe} DESTINATION lib/timing )
endforeach()
if( ${NUMBER_OF_CPU} LESS 4 )
    set( _p 1 )
    set( _q 1 )
elseif( ${NUMBER_OF_CPU} LESS 16 )
    set( _p 2 )
    set( _q 2 )
else()
    set( _p 4 )
    set( _q 4 )
endif()
math( EXPR _pq "${_p} * ${_q}" )

# Add basic tests. To run them type: 'make test' in build directory
add_test( NAME "timing_pdgesvd:A" COMMAND mpirun -np ${_pq} timing_pdgesvd --nprow ${_p} --npcol ${_q} --nb 64 --n_range 512:512:512 --niter 3 --check --qwmr )
add_test( NAME "timing_pdgesvd:B" COMMAND mpirun -np ${_pq} timing_pdgesvd --nprow ${_p} --npcol ${_q} --nb 64 --n_range 512:512:512 --niter 3 --check --qwdc )
add_test( NAME "timing_pdgesvd:C" COMMAND mpirun -np ${_pq} timing_pdgesvd --nprow ${_p} --npcol ${_q} --nb 64 --n_range 512:512:512 --niter 3 --check --qwel )
add_test( NAME "timing_pdgeqdwh:A" COMMAND mpirun -np ${_pq} timing_pdgeqdwh --nprow ${_p} --npcol ${_q} --nb 64 --n_range 512:512:512 --niter 3 --check )
