#!/usr/bin/env python3
"""Simple coverage analysis for crop/api.py"""

import re


def analyze_coverage():
    # Read the API file
    with open("src/crop/api.py", "r") as f:
        api_content = f.read()

    # Read the test file
    with open("tests/test_crop.py", "r") as f:
        test_content = f.read()

    # Find all function definitions in api.py
    function_pattern = r"^def\s+(\w+)\s*\("
    functions = re.findall(function_pattern, api_content, re.MULTILINE)

    print("=== CROP API Coverage Analysis ===\n")

    tested_functions = []
    untested_functions = []

    # Check which functions have tests
    for func in functions:
        # Look for test functions that reference this function
        test_patterns = [
            f"test.*{func}",  # test function names
            f"from.*import.*{func}",  # direct imports
            f"{func}\\(",  # function calls
        ]

        has_test = any(re.search(pattern, test_content, re.IGNORECASE) for pattern in test_patterns)

        if has_test:
            tested_functions.append(func)
        else:
            untested_functions.append(func)

    # Calculate coverage
    total_functions = len(functions)
    tested_count = len(tested_functions)
    coverage_percent = (tested_count / total_functions) * 100 if total_functions > 0 else 0

    print(f"Total Functions: {total_functions}")
    print(f"Tested Functions: {tested_count}")
    print(f"Untested Functions: {len(untested_functions)}")
    print(f"Function Coverage: {coverage_percent:.1f}%\n")

    print("TESTED FUNCTIONS:")
    for func in sorted(tested_functions):
        print(f"  ✓ {func}")

    print("\nUNTESTED FUNCTIONS:")
    for func in sorted(untested_functions):
        print(f"  ✗ {func}")

    # Analyze lines of code
    lines = api_content.split("\n")
    total_lines = len(lines)

    # Count non-empty, non-comment lines
    code_lines = [line for line in lines if line.strip() and not line.strip().startswith("#")]
    total_code_lines = len(code_lines)

    print("\n=== Line Analysis ===")
    print(f"Total Lines: {total_lines}")
    print(f"Code Lines (non-empty, non-comment): {total_code_lines}")

    # Estimate untested lines (very rough approximation)
    # Main untested function is run_crop_algorithm which is quite large
    if "run_crop_algorithm" in untested_functions:
        print("\nNOTE: run_crop_algorithm is a large function (~80+ lines)")
        print("Estimated untested code lines: ~80-100 lines")
        print(f"Estimated line coverage: ~{((total_code_lines - 90) / total_code_lines) * 100:.1f}%")


if __name__ == "__main__":
    analyze_coverage()
