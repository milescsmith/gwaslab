import os
import sys
import unittest

import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab.info.g_Log import Log
from gwaslab.qc.qc_fix_sumstats import flip_by_inverse


class TestFlipByInverse(unittest.TestCase):
    def test_flip_by_inverse_uses_original_ci_bounds_for_or_and_hr(self):
        sumstats = pd.DataFrame(
            {
                "OR": [1.10, 0.90, 1.25],
                "OR_95L": [0.80, 0.70, 1.05],
                "OR_95U": [1.40, 1.20, 1.60],
                "HR": [1.05, 1.20, 0.95],
                "HR_95L": [0.90, 1.00, 0.80],
                "HR_95U": [1.20, 1.50, 1.10],
            }
        )
        matched_index = pd.Series([True, False, True], index=sumstats.index)

        original_or_95l = sumstats["OR_95L"].copy()
        original_or_95u = sumstats["OR_95U"].copy()
        original_hr_95l = sumstats["HR_95L"].copy()
        original_hr_95u = sumstats["HR_95U"].copy()

        result = flip_by_inverse(
            sumstats=sumstats.copy(),
            matched_index=matched_index,
            log=Log(),
            verbose=False,
            factor=1,
        )

        expected_or_95u = original_or_95u.copy()
        expected_or_95l = original_or_95l.copy()
        expected_hr_95u = original_hr_95u.copy()
        expected_hr_95l = original_hr_95l.copy()

        expected_or_95u.loc[matched_index] = 1 / original_or_95l.loc[matched_index]
        expected_or_95l.loc[matched_index] = 1 / original_or_95u.loc[matched_index]
        expected_hr_95u.loc[matched_index] = 1 / original_hr_95l.loc[matched_index]
        expected_hr_95l.loc[matched_index] = 1 / original_hr_95u.loc[matched_index]

        np.testing.assert_allclose(result["OR_95U"], expected_or_95u)
        np.testing.assert_allclose(result["OR_95L"], expected_or_95l)
        np.testing.assert_allclose(result["HR_95U"], expected_hr_95u)
        np.testing.assert_allclose(result["HR_95L"], expected_hr_95l)


if __name__ == "__main__":
    unittest.main()
