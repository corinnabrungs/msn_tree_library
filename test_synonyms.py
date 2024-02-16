import ast
from unittest import TestCase
import pandas_utils as pu
from synonyms import parse_synonyms_to_list


class Test(TestCase):
    def test_get_all_synonyms(self):
        semicolon = "name;name2;hello"
        likelist = '["name", "name2", "hello"]'

        result = ["name", "name2", "hello"]

        self.assertEqual(result, parse_synonyms_to_list(semicolon))
        self.assertEqual([], parse_synonyms_to_list(None))
        self.assertEqual([], parse_synonyms_to_list(""))
        self.assertEqual(result, parse_synonyms_to_list(likelist))
