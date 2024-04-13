include_directories         (${CMAKE_SOURCE_DIR}/src)

add_executable              (trimmer                       trimmer.cpp)

target_link_libraries       (trimmer                       ${libraries})

set_property                (TARGET trimmer                PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)

install                     (TARGETS trimmer
                            DESTINATION ${CMAKE_INSTALL_PREFIX}/src/trimmer/
                            PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                            GROUP_READ GROUP_WRITE GROUP_EXECUTE
                            WORLD_READ WORLD_WRITE WORLD_EXECUTE)
