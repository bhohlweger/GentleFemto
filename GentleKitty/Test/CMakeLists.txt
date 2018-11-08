find_package(Boost COMPONENTS unit_test_framework REQUIRED)
include_directories (${Boost_INCLUDE_DIRS})

add_executable(lambdaParam testLambdaParam.cxx)
target_link_libraries(lambdaParam ${Boost_LIBRARIES} tidyCATS)
add_test(NAME lambdaParam COMMAND lambdaParam)