"""Ensure the installed wheel's leafblower package is found, not the local
source tree which lacks the compiled _leafblower extension."""
import sys
import pathlib

_here = str(pathlib.Path(__file__).parent)
if _here in sys.path:
    sys.path.remove(_here)
