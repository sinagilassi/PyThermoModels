from setuptools import setup, find_packages
import codecs
import os

here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as fh:
    long_description = "\n" + fh.read()

APP_NAME = 'PyThermoModels'
VERSION = '1.2.0'
AUTHOR = 'Sina Gilassi'
EMAIL = "<sina.gilassi@gmail.com>"
DESCRIPTION = 'A Python package designed for the calculation of thermodynamic properties using various well-established models.'
LONG_DESCRIPTION = "PyThermoModels is an open-source Python package designed to facilitate thermodynamic modeling and calculations. This package provides a comprehensive and user-friendly interface to popular thermodynamic models, enabling quick and accurate estimation of key properties."

# Setting up
setup(
    name=APP_NAME,
    version=VERSION,
    author=AUTHOR,
    author_email=EMAIL,
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=long_description,
    packages=find_packages(exclude=['tests', '*.tests', '*.tests.*']),
    include_package_data=True,  # Make sure to include non-Python files
    # Add both config and data files
    package_data={'': ['config/*.yml', 'config/*.csv', 'data/*.csv',
                       'data/*.yml', 'plugin/*.yml', 'plugin/*.csv']},
    license='MIT',
    license_files=[],
    install_requires=['pandas', 'numpy',
                      'PyYAML', 'PyCUC', 'scipy',
                      'PyThermoDB', 'PyThermoLinkDB'],
    extras_require={
        "plotting": ["matplotlib"],
    },
    keywords=[
        'python', 'chemical-engineering', 'equation-of-state',
        'thermodynamic-properties', 'activity-coefficient',
        'thermodynamic-models', 'NRTL', 'UNIQUAC', 'SRK', 'PR'
    ],
    classifiers=[
        "Development Status :: Planning",
        "Intended Audience :: Education",
        "Topic :: Scientific/Engineering",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ],
    python_requires='>=3.10',
)
