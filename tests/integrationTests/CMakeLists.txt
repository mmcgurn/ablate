add_executable(integrationTests "")
target_link_libraries(integrationTests PUBLIC gtest ablateLibrary testingResources)

target_sources(integrationTests
        PRIVATE
        tests.cpp
        main.cpp
        )

gtest_discover_tests(integrationTests
        # set a working directory so your project root so that you can find test data via paths relative to the project root
        WORKING_DIRECTORY ${PROJECT_DIR}
        PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_DIR}"
        )

add_custom_command(
        TARGET integrationTests
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/outputs/
        $<TARGET_FILE_DIR:integrationTests>/outputs
)

add_custom_command(
        TARGET integrationTests
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/inputs/
        $<TARGET_FILE_DIR:integrationTests>/inputs
)

