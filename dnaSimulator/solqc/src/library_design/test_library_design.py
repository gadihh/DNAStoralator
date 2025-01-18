from __future__ import print_function
import unittest
import pandas as pd

from src.library_design.library_design import LibraryDesign


def get_library_design_object():
    design_df = pd.read_csv("/Users/yoavorlev/Desktop/library/data/Twist_zohar_design.csv", nrows=5)
    return LibraryDesign(design_df)


class TestUM(unittest.TestCase):
    def setUp(self):
        pass

    def test_get_design_longest_read(self):
        lb_object = get_library_design_object()
        longest_read = lb_object.get_design_longest_sequence_len()
        print(longest_read)
        self.assertEqual(longest_read, 148)


if __name__ == '__main__':
    unittest.main()