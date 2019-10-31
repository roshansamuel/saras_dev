include_directories (
    ${PARENT_DIR}/src/solvers/hydro
    ${PARENT_DIR}/src/solvers/scalar
)

set (EXECUTABLE_OUTPUT_PATH ${PARENT_DIR})

set (HYDRO_SOURCE hydro/hydro.cc hydro/hydro_d2.cc hydro/hydro_d3.cc)

set (SCALAR_SOURCE scalar/scalar.cc scalar/scalar_d2.cc scalar/scalar_d3.cc)

set (SOURCES main.cc ${HYDRO_SOURCE} ${SCALAR_SOURCE})

add_executable (saras ${SOURCES})

target_link_libraries(saras field grid parser probes reader writer boundary parallel poisson force yaml-cpp hdf5 debug /usr/local/lib/libblitz.a)
