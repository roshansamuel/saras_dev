include_directories (
    ${PARENT_DIR}/src/solvers/hydro
    ${PARENT_DIR}/src/tests/global
    ${PARENT_DIR}/src/tests/units
)

set (EXECUTABLE_OUTPUT_PATH ${PARENT_DIR})

set (HYDRO_SOURCE ${PARENT_DIR}/src/solvers/hydro/hydro.cc ${PARENT_DIR}/src/solvers/hydro/hydro2.cc ${PARENT_DIR}/src/solvers/hydro/hydro3.cc)

set (TEST_SOURCE units/field_test.cc units/differ_test.cc units/nlin_test.cc units/poisson_test.cc units/hydro_test.cc)

set (UNIT_SOURCE global/unittest.cc)

set (SOURCES main.cc ${UNIT_SOURCE} ${HYDRO_SOURCE} ${TEST_SOURCE})

add_executable(saras_test ${SOURCES})

target_link_libraries(saras_test field grid parser probes reader writer parallel poisson yaml-cpp hdf5 debug /usr/local/lib/libblitz.a)
