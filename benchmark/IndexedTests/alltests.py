#!/usr/bin/python
# vim: set expandtab ts=4
import unittest
from concurrencytest import ConcurrentTestSuite, fork_for_tests
import IndexedTest
import ChangeOfBaseTests

if  __name__ == '__main__':
    s= unittest.TestLoader().loadTestsFromTestCase(IndexedTest.IndexedTest)
    suites=[ChangeOfBaseTests.suite()]
    for suite in suites:
        s.addTests(suite)
    # Run same tests across 16 processes
    concurrent_suite = ConcurrentTestSuite(s, fork_for_tests(16))
    runner = unittest.TextTestRunner()
    runner.run(concurrent_suite)
