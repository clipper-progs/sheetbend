cmake_minimum_required(VERSION 3.1)
project(SHEETBEND C CXX)

set(CMAKE_CXX_STANDARD 11)

find_package(CCP4 REQUIRED mmdb2 ccp4c clipper-core clipper-mmdb clipper-minimol
                           clipper-contrib clipper-ccp4 clipper-cif)

include_directories(${CCP4_INCLUDE_DIRS}
                    ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(csheetbend
               csheetbend.cpp)

install(TARGETS csheetbend
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION "${LIB_INSTALL_DIR}"
        LIBRARY DESTINATION "${LIB_INSTALL_DIR}")
install(TARGETS csheetbend DESTINATION bin)

if (UNIX)
  install(CODE "
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink csheetbend
            \"\$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/csheetbend\")
        ")
endif()
