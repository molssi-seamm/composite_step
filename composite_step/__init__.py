# -*- coding: utf-8 -*-

"""
composite_step
A SEAMM plug-in for composite methods: energies from recipes of ordinary calculations
"""

# Bring up the classes so that they appear to be directly in
# the composite_step package.

from .composite import Composite  # noqa: F401, E501
from .composite_parameters import CompositeParameters  # noqa: F401, E501
from .composite_step import CompositeStep  # noqa: F401, E501
from .tk_composite import TkComposite  # noqa: F401, E501

from .metadata import metadata  # noqa: F401

# Handle versioneer
from ._version import get_versions

__author__ = "Paul Saxe"
__email__ = "psaxe@molssi.org"
versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
