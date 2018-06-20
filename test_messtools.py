"""This module contains tests for messtools python bindings
"""
import numpy
import messtools


def test__mass():
    """test messtools.mass
    """
    assert numpy.allclose(messtools.mass('C'), 21874.661832)
