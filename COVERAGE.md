# Code Coverage for CROP

This project now includes comprehensive code coverage analysis using `pytest-cov` and `coverage`.

## Quick Start

### Run Coverage Analysis

```bash
# Simple coverage run (API utility functions only)
uv run python -m pytest tests/test_crop.py::TestApiUtilityFunctions --cov=src --cov-report=term-missing

# Full coverage analysis (excluding solver-dependent tests)
uv run python -m pytest tests/ --cov=src --cov-report=term-missing -k "not test_crop_algorithm_integration and not test_complete_crop_workflow"

# Use the comprehensive coverage script
python run_coverage.py
```

### View Coverage Reports

After running coverage, you can view the results in multiple formats:

1. **Terminal Output**: Shows coverage percentages and missing lines
2. **HTML Report**: Open `htmlcov/index.html` in your browser for detailed interactive coverage
3. **XML Report**: `coverage.xml` for CI/CD integration

## Coverage Configuration

The coverage settings are configured in `pyproject.toml`:

- **Source**: `src/` directory
- **Branch Coverage**: Enabled for comprehensive analysis
- **Reports**: Terminal, HTML, and XML formats
- **Exclusions**: Test files, setup files, and common patterns

## Current Coverage

As of the latest run:

### API Utility Functions (tests/test_crop.py::TestApiUtilityFunctions)
- **25 tests** covering all utility functions in `src/crop/api.py`
- **100% function coverage** for API utility functions
- **Comprehensive edge case testing** including error handling

### Overall Project Coverage
- **src/crop/api.py**: ~62% coverage (main focus)
- **src/crop/__init__.py**: 100% coverage
- **Total**: ~53-58% overall coverage

## What's Covered

✅ **Fully Tested**:
- `build_phenotype_conditions` - 8 comprehensive test cases
- `get_growth_conditions` - 3 test scenarios
- `get_nogrowth_conditions` - 3 test scenarios  
- `get_carbon_source` - Unit tested with mocks
- All bound calculation functions (7 functions)
- Error handling and edge cases
- Data transformation and validation

⚠️ **Partially Covered**:
- `run_crop_algorithm` - Integration tested but requires specific solvers
- Complex CVXPY constraint generation functions

❌ **Not Covered**:
- CLI functionality (`src/crop/cli.py`)
- Version module dynamic parts
- Main entry point

## Adding New Tests

When adding new functionality, ensure you:

1. **Write unit tests** for all new functions
2. **Include edge cases** and error handling
3. **Use mocks** for complex dependencies
4. **Run coverage** to verify completeness

Example test structure:
```python
def test_new_function():
    # Test normal case
    result = new_function(valid_input)
    assert result == expected_output
    
    # Test edge case
    with pytest.raises(ValueError):
        new_function(invalid_input)
```

## CI/CD Integration

The `coverage.xml` file can be used with CI/CD systems:

```yaml
# Example GitHub Actions
- name: Run tests with coverage
  run: |
    uv run python -m pytest --cov=src --cov-report=xml
    
- name: Upload coverage to Codecov
  uses: codecov/codecov-action@v3
  with:
    file: ./coverage.xml
```

## Tools Used

- **pytest**: Test framework
- **pytest-cov**: Coverage plugin for pytest
- **coverage**: Core coverage measurement tool
- **uv**: Package and environment management

## Troubleshooting

### Solver Errors
Some integration tests require specific solvers (like GUROBI). Use the filtering options to exclude these:
```bash
-k "not test_crop_algorithm_integration and not test_complete_crop_workflow"
```

### Missing Coverage
If coverage seems low, check:
1. Are all test modules being discovered?
2. Are imports working correctly?
3. Are there unused/dead code sections?

### HTML Report Not Generated
Ensure you have write permissions in the project directory and that the `htmlcov` directory can be created.