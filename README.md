# chisq_test_wrapper
A friendly, automated chi-square test function which takes care of post-hoc tests and multiple comparisons.

This module contains functions that wrap the standard chi2_contingency test from scipy.stats.,
as well as additional corrections for multiple comparisons and post-hoc tests.

The usage is as simple as chisq_and_posthoc_corrected(df).
