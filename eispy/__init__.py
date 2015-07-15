# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This is an Astropy affiliated package.
"""

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from __future__ import absolute_import
# ----------------------------------------------------------------------------

try:
    from .version import version as __version__
except ImportError:
    __version__ = ''
    
try:
    from .version import githash as __githash__
except ImportError:
    __githash__ = ''