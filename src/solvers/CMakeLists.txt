include_directories (
    ${PARENT_DIR}/src/solvers/hydro
    ${PARENT_DIR}/src/solvers/thermal
)

set (EXECUTABLE_OUTPUT_PATH ${PARENT_DIR})

set (HYDRO_SOURCE hydro/hydro.cc hydro/hydro2.cc hydro/hydro3.cc)

set (THERMAL_SOURCE thermal/thermal.cc thermal/thermal2.cc thermal/thermal3.cc)

set (SOURCES main.cc ${HYDRO_SOURCE} ${THERMAL_SOURCE})

add_executable (saras ${SOURCES})

target_link_libraries(saras field grid parser probes reader writer boundary parallel poisson yaml-cpp hdf5 debug /usr/local/lib/libblitz.a)
