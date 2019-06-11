#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

import functools
import nose

def expected_failure(test):
    @functools.wraps(test)
    def inner(*args, **kwargs):
        try:
            test(*args, **kwargs)
        except Exception:
            raise nose.SkipTest
        else:
            raise AssertionError('Failure expected')
    return inner
