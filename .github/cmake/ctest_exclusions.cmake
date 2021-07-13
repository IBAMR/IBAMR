# List tests here whose results should be omitted.
set(test_exclusions
)

# Platform specific exclusions:
if ("$ENV{CMAKE_CONFIGURATION}" MATCHES "fedora")
  list(APPEND test_exclusions
    # Comment on why the test fails.
    # Regular expression matching test: "^RenderMesh$"
  )
endif ()

string(REPLACE ";" "|" test_exclusions "${test_exclusions}")
if (test_exclusions)
  set(test_exclusions "(${test_exclusions})")
endif ()
