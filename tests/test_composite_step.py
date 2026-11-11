#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `composite_step` package."""

import pytest  # noqa: F401
import composite_step  # noqa: F401


def test_construction():
    """Just create an object and test its type."""
    result = composite_step.Composite()
    assert str(type(result)) == "<class 'composite_step.composite.Composite'>"
