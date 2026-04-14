# only a function test, not correctness :-|
function(test_fmesher file)
    add_test(NAME fmesher_${file}.setup
        COMMAND "${CMAKE_COMMAND}"
        -E copy "${CMAKE_CURRENT_LIST_DIR}/${file}" .
        )
    add_test(NAME fmesher_${file}
        COMMAND fmesher-bin "${file}"
        )
    set_tests_properties(fmesher_${file} PROPERTIES
        DEPENDS fmesher_${file}.setup
        LABELS "magnetics;mesher"
        )
endfunction()

function(test_fmesher_virtual_gap file)
    add_test(NAME fmesher_${file}.setup
        COMMAND "${CMAKE_COMMAND}"
        -E copy "${CMAKE_CURRENT_LIST_DIR}/${file}" .
        )
    add_test(NAME fmesher_${file}
        COMMAND fmesher-bin "${file}"
        )
    set_tests_properties(fmesher_${file} PROPERTIES
        DEPENDS fmesher_${file}.setup
        LABELS "magnetics;mesher;virtual-gap"
        )
    add_test(NAME fmesher_${file}.check
        COMMAND "${CMAKE_COMMAND}"
        -D FEM_FILE=${file}
        -P "${CMAKE_CURRENT_LIST_DIR}/check_virtual_gap_rewrite.cmake"
        )
    set_tests_properties(fmesher_${file}.check PROPERTIES
        DEPENDS fmesher_${file}
        LABELS "magnetics;mesher;virtual-gap"
        )
endfunction()

function(test_fmesher_virtual_gap_failure file)
    add_test(NAME fmesher_${file}.setup
        COMMAND "${CMAKE_COMMAND}"
        -E copy "${CMAKE_CURRENT_LIST_DIR}/${file}" .
        )
    add_test(NAME fmesher_${file}
        COMMAND fmesher-bin "${file}"
        )
    set_tests_properties(fmesher_${file} PROPERTIES
        DEPENDS fmesher_${file}.setup
        WILL_FAIL TRUE
        LABELS "magnetics;mesher;virtual-gap"
        )
endfunction()

test_fmesher(Temp.fem)
test_fmesher(split_seg_err_test.fem)
test_fmesher_virtual_gap(virtual_gap_line_simple.fem)
test_fmesher_virtual_gap_failure(virtual_gap_periodic_guard.fem)
# vi:expandtab:tabstop=4 shiftwidth=4:
