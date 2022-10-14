# Required declarations
cmake_minimum_required(VERSION 3.12.0)
project(minimum_disc CXX)
set(CMAKE_CXX_VERSION 20)

# Collect sources and includes
include_directories(inc)
file(GLOB_RECURSE sources CONFIGURE_DEPENDS "src/*.cpp")

# Add external dependencies.
# Note: find_package with REQUIRED will probably not work with emscripten,
# so we'll need to come up with something else => FetchContent for example.
find_package(Eigen3)
if(NOT ${Eigen3_FOUND})
    include(FetchContent)
    FetchContent_Declare(eigen3 GIT_REPOSITORY "https://gitlab.com/libeigen/eigen" GIT_TAG 3.4.0)
    FetchContent_MakeAvailable(eigen3)

    set(EIGEN3_INCLUDE_DIR "${eigen3_SOURCE_DIR}")
endif()

# Output the include structure to help the language server,
# and symlink it to the source directory.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink
                "${CMAKE_BINARY_DIR}/compile_commands.json"
                "${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json")

# Build binaries
add_library(${PROJECT_NAME} STATIC ${sources})
target_include_directories(${PROJECT_NAME} PUBLIC inc "${EIGEN3_INCLUDE_DIR}")

# Build tests
file(GLOB_RECURSE testSources test/*.cpp)
add_executable(${PROJECT_NAME}_testrunner ${testSources})
target_link_libraries(${PROJECT_NAME}_testrunner ${PROJECT_NAME})
