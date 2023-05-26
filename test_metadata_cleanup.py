from unittest import TestCase
import pandas as pd

import metadata_cleanup


def check_smilesdf(df):
    clean = metadata_cleanup.ensure_smiles_column(df)
    assert clean is not None
    assert "smiles" in clean.columns
    assert clean.at[0, "smiles"] is not None
    assert clean.at[0, "smiles"] == "iso"
    assert clean.at[2, "smiles"] == "can"
    assert clean.at[3, "smiles"] == "iso"


class Test(TestCase):

    def test_ensure_smiles_column(self):
        check_smilesdf(
            pd.DataFrame(
                {
                    "smiles": [None, "lfjlajdfi", None, ""],
                    "isomeric_smiles": ["iso", None, None, "iso"],
                    "canonical_smiles": ["can", "can", "can", "can"]
                }
            )
        )

        check_smilesdf(
            df=pd.DataFrame(
                {
                    "isomeric_smiles": ["iso", None, None, "iso"],
                    "canonical_smiles": ["can", "can", "can", "can"]
                }
            )
        )

        clean = metadata_cleanup.ensure_smiles_column(pd.DataFrame(
            {
                "inchi": ["", None]
            }
        ))
        assert clean is not None
        assert "smiles" in clean
