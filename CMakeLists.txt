cmake_minimum_required(VERSION 3.21)

project(FactorioAssignment LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 1. Include the FetchContent module to download dependencies
include(FetchContent)

# 2. Declare the nlohmann-json dependency
FetchContent_Declare(
  nlohmann_json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG        v3.11.3 # Use a specific version for determinism
)
# Make it available (this downloads it)
FetchContent_MakeAvailable(nlohmann_json)


# 3. Declare the Google OR-Tools dependency
FetchContent_Declare(
  ortools
  GIT_REPOSITORY https://github.com/google/or-tools.git
  GIT_TAG        v9.9    # Use the same version for determinism
)
# IMPORTANT: Tell ortools to download its own pre-built binaries
# instead of trying to build its dependencies (like Protobuf) from source.
# This is MUCH faster and avoids errors.
set(BUILD_DEPS ON CACHE BOOL "Build OR-Tools dependencies")

# Make it available (this downloads it and configures it)
FetchContent_MakeAvailable(ortools)


# 4. Add your part2_assignment subdirectory
add_subdirectory(part2_assignment)