# Define source files
set( POLAR_SRC pdgenm2.c
               mellipj.c
               mellipke.c
               pdgezolopd.c
               pdgeqdwh.c
               )


list( APPEND CMAKE_SHARED_LINKER_FLAGS "-Wl,--enable-new-dtags")
# Transform the list into a space separate string
string(REPLACE ";" " " CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}")

# Build POLAR
add_library( ${POLAR_LIBNAME} ${POLAR_SRC} )
# Libraries to link
target_link_libraries( ${POLAR_LIBNAME} ${MPI_C_LIBRARIES} ${SCALAPACK_LIBRARIES} ${EXTRA_LIBS} )

# Installation
install( TARGETS ${POLAR_LIBNAME} DESTINATION lib )
