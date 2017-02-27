# Define source files
set( QDWH_SRC pdgenm2.c
              pdgeqdwh.c
              pdgeqsvd.c )


list( APPEND CMAKE_SHARED_LINKER_FLAGS "-Wl,--enable-new-dtags")
# Transform the list into a space separate string
string(REPLACE ";" " " CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}")

# Build QDWH
add_library( ${QDWH_LIBNAME} ${QDWH_SRC} )
# Libraries to link
target_link_libraries( ${QDWH_LIBNAME} ${MPI_C_LIBRARIES} ${SCALAPACK_LIBRARIES} ${ELPA_LIBRARIES} ${EXTRA_LIBS} )

# Add ELPA depencency, to build ELPA first
if( NOT ELPA_FOUND )
    add_dependencies(${QDWH_LIBNAME} ELPA)
endif()

# Installation
install( TARGETS ${QDWH_LIBNAME} DESTINATION lib )
