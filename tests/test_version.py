import propka
import re

def test_version_exists():
    assert hasattr(propka, '__version__')

def test_version():
    assert (re.match(r'[0-9]+\.[0-9]+', propka.__version__) or
            propka.__version__.startswith('0+untagged'))

