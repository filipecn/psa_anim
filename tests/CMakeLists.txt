add_executable(psa_anim_tests
        main.cpp
        # dsl_data_tests.cpp
        # stl_foam_test.cpp
        interpolation_tests.cpp
        )
target_include_directories(psa_anim_tests PUBLIC
        ${CATCH2_INCLUDES}
        ${PSA_ANIM_INCLUDES}
        )
        add_definitions(-DASSETS_PATH="${CMAKE_CURRENT_LIST_DIR}/assets")

target_link_options(psa_anim_tests INTERFACE --coverage)

target_link_libraries(psa_anim_tests PRIVATE
        Catch2::Catch2WithMain
        Eigen3::Eigen
        gcov
        ${PSA_ANIM_LIBRARIES}
        )

add_custom_target(test_psa_anim
        COMMAND lcov --directory . --zerocounters
        COMMAND psa_anim_tests DEPENDS psa_anim_tests
        COMMAND lcov --directory . --capture -output-file cov.info
        COMMAND lcov --remove cov.info '/usr/*' --output-file cov.info
        COMMAND lcov --remove cov.info '*/ext/*' --output-file cov.info
        COMMAND lcov --remove cov.info '*/tests/*' --output-file cov.info
        COMMAND genhtml cov.info --output-directory out
        )
