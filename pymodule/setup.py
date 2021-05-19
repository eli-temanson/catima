from pathlib import Path

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

DIR = Path(__file__).parents[0]

SRC = [str(DIR/'pycatima.cpp')]+[str(fname.resolve()) for fname in DIR.glob('../*.cpp')]
SRCG = [str(fname.resolve()) for fname in DIR.glob('../global/*.c')]
SRC += SRCG
print (SRC)
example_module = Pybind11Extension(
    'pycatima',
    SRC,
    include_dirs=['../build/include','../global']
)

setup(
    name='pycatima',
    version=1.54,
    author='Andrej Prochazka',
    author_email='hrocho@vodacionline.sk',
    description='python interface to catima library',
    ext_modules=[example_module],
    cmdclass={"build_ext": build_ext},
)
