add_executable(belts main.cpp)

target_link_libraries(belts
  PRIVATE
    ortools::ortools
    nlohmann_json::nlohmann_json
)

set_target_properties(belts PROPERTIES
    CXX_STANDARD ${CMAKE_CXX_STANDARD}
    CXX_STANDARD_REQUIRED ${CMAKE_CXX_STANDARD_REQUIRED}
)