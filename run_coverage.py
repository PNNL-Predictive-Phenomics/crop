#!/usr/bin/env python3
"""
Code coverage script for the CROP project.

This script runs the test suite with code coverage analysis and generates
comprehensive reports including HTML, XML, and terminal output.
"""

import subprocess
import sys
from pathlib import Path


def run_command(cmd, description):
    """Run a command and handle errors."""
    print(f"\n🚀 {description}")
    print(f"Running: {' '.join(cmd)}")
    print("-" * 60)
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=False)
        print(f"✅ {description} completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ {description} failed with exit code {e.returncode}")
        return False


def main():
    """Run comprehensive code coverage analysis."""
    print("🧪 CROP Code Coverage Analysis")
    print("=" * 60)
    
    # Commands to run
    commands = [
        # Run API utility function tests with coverage
        (
            ["uv", "run", "python", "-m", "pytest", 
             "tests/test_crop.py::TestApiUtilityFunctions", 
             "--cov=src", 
             "--cov-report=term-missing", 
             "--cov-report=html", 
             "--cov-report=xml",
             "--cov-branch",
             "-v"],
            "API Utility Functions Coverage"
        ),
        
        # Run tests excluding problematic integration tests
        (
            ["uv", "run", "python", "-m", "pytest", 
             "tests/",
             "--cov=src", 
             "--cov-report=term-missing", 
             "--cov-report=html", 
             "--cov-report=xml",
             "--cov-branch",
             "-v",
             "-k", "not test_crop_algorithm_integration and not test_complete_crop_workflow"],
            "Full Test Suite Coverage (excluding solver-dependent tests)"
        ),
    ]
    
    # Run coverage commands
    success_count = 0
    for cmd, description in commands:
        if run_command(cmd, description):
            success_count += 1
    
    # Summary
    print("\n" + "=" * 60)
    print("📊 COVERAGE ANALYSIS SUMMARY")
    print("=" * 60)
    
    if success_count == len(commands):
        print("✅ All coverage analysis completed successfully!")
        
        # Check if HTML coverage report exists
        html_report = Path("htmlcov/index.html")
        if html_report.exists():
            print(f"\n📁 HTML Coverage Report: {html_report.absolute()}")
            print("   Open this file in your browser to view detailed coverage")
        
        # Check if XML coverage report exists
        xml_report = Path("coverage.xml")
        if xml_report.exists():
            print(f"\n📄 XML Coverage Report: {xml_report.absolute()}")
            print("   Use this file for CI/CD integration")
            
        print("\n🎯 Coverage Analysis Complete!")
        print("Check the terminal output above for detailed coverage metrics.")
        
    else:
        print(f"❌ {len(commands) - success_count} out of {len(commands)} coverage analyses failed")
        sys.exit(1)


if __name__ == "__main__":
    main()