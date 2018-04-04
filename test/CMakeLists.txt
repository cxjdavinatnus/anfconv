add_executable(testpolybori
    testpolybori.cpp
)

target_link_libraries(testpolybori
    libpolybori
    libcudd
)