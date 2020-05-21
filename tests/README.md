# Testing PROPKA

These tests assume that PROPKA is installed as a module on your system.
If you are running in a virtual environment and want to make changes to your
code, module installation accomplished by
```
pip install -e .
```
from the top level of the PROPKA source directory.

Once PROPKA is available as a module, the tests can be run by
```
python -m pytest tests
```
either in the top-level directory or `tests` subdirectory.
