add_executable(test_app
    test.cpp
    greedy_b_matching.cpp
    anstee_b_matching.cpp
)
target_link_libraries(test_app GreeBaLib ortools::ortools)
