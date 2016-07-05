"""
Copyright (C) 2014 Jakub Krajniak <jkrajniak@gmail.com>

This file is part of BondMatcher.

BondMatcher is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

BondMatcher - unit tests for the analysis part
"""

import unittest

class BaseTestCase(unittest.TestCase):

  def assertSameElements(self, one, two):
    """Checks if two iterable elements are the same. The order does not matter."""
    missing = []

    for el_one in one:
      if el_one not in two:
        missing.append(el_one)

    self.assertEqual(missing, [])

  def assertSameElementsNestedList(self, one, two):
    count = 0
    for el_one in one:
      for el_two in two:
        if set(el_one) == set(el_two):
          count += 1

    self.assertEqual(count, len(one))

class TestBaseTestCase(BaseTestCase):

  def testAssertSameElementsOk(self):
    one = [1, 2, 3, 4]
    two = [3, 4, 1, 2]

    self.assertSameElements(one, two)

  def testAsssertsameElementsFail(self):
    one = [1, 2, 3]
    two = [1, 2, 4]

    try:
      self.assertSameElements(one, two)
    except Exception as ex:
      self.assertIsInstance(ex, AssertionError)


if __name__ == '__main__':
  unittest.main()
