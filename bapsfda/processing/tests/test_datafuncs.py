import numpy as np
from bapsfda.processing import abs2


class TestProcessing:
    def test_abs2(self):
        a = np.array([1, 2, 3, 4])
        b = np.array([1 + 1j, 2 + 2j, 3 + 3j, 4 + 4j])
        assert np.array_equal(abs2(a), np.array([1, 4, 9, 16])), f"fail for real array"
        assert np.array_equal(
            abs2(b), np.array([2, 8, 18, 32])
        ), f"fail for complex array"
