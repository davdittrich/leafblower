"""Root conftest: remove python/ src tree from sys.path so the installed wheel
is imported instead of the local source directory (which lacks _leafblower.so)."""
import sys
import pathlib

_python_src = str(pathlib.Path(__file__).parent / "python")
if _python_src in sys.path:
    sys.path.remove(_python_src)
