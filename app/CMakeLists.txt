add_executable(test_app
    test.cpp
    greedy_b_matching.cpp
    anstee_b_matching.cpp
    reduction_b_matching.cpp
    milp_b_matching.cpp
)
target_link_libraries(test_app GreeBaLib ortools::ortools)

if(BUILD_TESTING)
    add_executable(reduction_tests
        reduction_tests.cpp
        reduction_b_matching.cpp
        anstee_b_matching.cpp
        milp_b_matching.cpp
    )
    target_link_libraries(reduction_tests GreeBaLib ortools::ortools)
    add_test(NAME reduction_b_matching COMMAND reduction_tests)
endif()
