add_executable(factory main.cpp)

target_link_libraries(factory
  PRIVATE
    ortools::ortools
    nlohmann_json::nlohmann_json
)

set_target_properties(factory PROPERTIES
    CXX_STANDARD ${CMAKE_CXX_STANDARD}
    CXX_STANDARD_REQUIRED ${CMAKE_CXX_STANDARD_REQUIRED}
)