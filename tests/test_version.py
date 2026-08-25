# -*- coding: utf-8 -*-

"""Comprehensive tests for version module."""

import unittest
from subprocess import CalledProcessError
from unittest.mock import patch

from crop.version import VERSION, get_git_hash, get_version


class TestVersion(unittest.TestCase):
    """Test version module functions."""

    def test_version_type(self):
        """Test the version is a string."""
        version = get_version()
        self.assertIsInstance(version, str)

    def test_version_string_format(self):
        """Test version string has expected format."""
        version = get_version()
        # Should be at least the VERSION constant
        self.assertIn("0.0.1", version)

    def test_get_version_without_git_hash(self):
        """Test get_version returns version without git hash."""
        version = get_version(with_git_hash=False)
        self.assertIsInstance(version, str)
        self.assertEqual(version, VERSION)

    def test_get_version_with_git_hash_false(self):
        """Test explicit with_git_hash=False parameter."""
        version_without = get_version(with_git_hash=False)
        version_default = get_version()
        self.assertEqual(version_without, version_default)

    @patch("crop.version.check_output")
    def test_get_git_hash_success(self, mock_check_output):
        """Test get_git_hash returns hash on success."""
        mock_check_output.return_value = b"abc1234567890def\n"
        git_hash = get_git_hash()
        self.assertIsInstance(git_hash, str)
        self.assertEqual(git_hash, "abc12345")  # First 8 chars

    @patch("crop.version.check_output")
    def test_get_git_hash_strips_whitespace(self, mock_check_output):
        """Test get_git_hash properly strips whitespace and limits to 8 chars."""
        mock_check_output.return_value = b"  fedcba9876543210  \n"
        git_hash = get_git_hash()
        self.assertEqual(git_hash, "fedcba98")

    @patch("crop.version.check_output")
    def test_get_git_hash_error_handling(self, mock_check_output):
        """Test get_git_hash returns UNHASHED on error."""
        mock_check_output.side_effect = CalledProcessError(1, "git")
        git_hash = get_git_hash()
        self.assertEqual(git_hash, "UNHASHED")

    @patch("crop.version.get_git_hash")
    def test_get_version_with_git_hash(self, mock_get_git_hash):
        """Test get_version includes git hash when requested."""
        mock_get_git_hash.return_value = "abc12345"
        version = get_version(with_git_hash=True)
        self.assertIsInstance(version, str)
        self.assertIn("-abc12345", version)
        self.assertTrue(version.startswith(VERSION))

    @patch("crop.version.get_git_hash")
    def test_get_version_with_git_hash_unhashed(self, mock_get_git_hash):
        """Test get_version with UNHASHED git state."""
        mock_get_git_hash.return_value = "UNHASHED"
        version = get_version(with_git_hash=True)
        self.assertIn("-UNHASHED", version)

    def test_get_git_hash_returns_string(self):
        """Test get_git_hash always returns a string."""
        result = get_git_hash()
        self.assertIsInstance(result, str)

    def test_version_constant_exists(self):
        """Test that VERSION constant is defined."""
        self.assertIsInstance(VERSION, str)
        self.assertEqual(VERSION, "0.0.1-dev")
