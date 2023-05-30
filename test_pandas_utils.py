from unittest import TestCase

import pandas as pd

import pandas_utils


class Test(TestCase):
    def test_remove_empy_strings(self):
        df = pd.DataFrame(
            {
                "col": [None, "", "sads"]
            }
        )
        df = pandas_utils.remove_empty_strings(df, columns="col")
        assert df is not None
        assert len(df[df["col"].isnull()]) == 2
        assert df.at[1, "col"] is None
