"""Ensure the installed wheel's leafblower package is found, not the local
source tree which lacks the compiled _leafblower extension."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import sys
import pathlib

_here = str(pathlib.Path(__file__).parent)
if _here in sys.path:
    sys.path.remove(_here)
