include_directories         (${CMAKE_SOURCE_DIR}/src)

add_executable              (eval                       eval.cpp)

target_link_libraries       (eval                       ${libraries})

set_property                (TARGET eval                PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)

install                     (TARGETS eval
                            DESTINATION ${CMAKE_INSTALL_PREFIX}/src/eval/
                            PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                            GROUP_READ GROUP_WRITE GROUP_EXECUTE
                            WORLD_READ WORLD_WRITE WORLD_EXECUTE)
