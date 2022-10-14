# Required declarations
cmake_minimum_required(VERSION 3.12.0)
project(wasm_demo CXX)

set(CMAKE_CXX_STANDARD 20)

# Collect sources and includes
file(GLOB_RECURSE sources CONFIGURE_DEPENDS "src/*.cpp")

# Output the include structure to help the language server,
# and symlink it to the source directory.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink
                "${CMAKE_BINARY_DIR}/compile_commands.json"
                "${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json")

# Add included dependencies
add_subdirectory(external)

# Options
set(WASM_DEMO_COMPILER_ARGUMENTS "" CACHE STRING "Compiler arguments passed to the WASM wrappers")
set(WASM_DEMO_LINKER_ARGUMENTS "" CACHE STRING "Compiler arguments passed to the WASM wrappers")

# Build binaries
include_directories(inc)
add_library(${PROJECT_NAME} STATIC ${sources})
target_link_libraries(${PROJECT_NAME} PUBLIC minimum_disc)

# Copy HTML and JS sources
file(GLOB_RECURSE htmls src/*.html)
file(GLOB_RECURSE javascripts src/*.js)

set(webSources "")
list(APPEND webSources ${htmls})
list(APPEND webSources ${javascripts})

foreach(webSource ${webSources})
    get_filename_component(webSourceName "${webSource}" NAME)
    #configure_file(${html} "${CMAKE_BINARY_DIR}/${html_name}" COPYONLY)
    EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink
                    "${webSource}"
                    "${CMAKE_BINARY_DIR}/${webSourceName}")
endforeach()
