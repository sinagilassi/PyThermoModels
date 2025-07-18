from setuptools import setup, find_packages
import codecs
import os

here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as fh:
    long_description = "\n" + fh.read()

APP_NAME = 'PyThermoModels'
VERSION = '1.4.8'
AUTHOR = 'Sina Gilassi'
EMAIL = "<sina.gilassi@gmail.com>"
DESCRIPTION = 'A Python package designed for the calculation of thermodynamic properties using various well-established models.'
LONG_DESCRIPTION = "PyThermoModels is an open-source Python package designed to facilitate thermodynamic modeling and calculations. This package provides a comprehensive and user-friendly interface to popular thermodynamic models, enabling quick and accurate estimation of key properties."
HOME_PAGE = 'https://github.com/sinagilassi/PyThermoModels'
DOCUMENTATION = 'https://pythermomodels.readthedocs.io/en/latest/'

# Setting up
setup(
    name=APP_NAME,
    version=VERSION,
    author=AUTHOR,
    author_email=EMAIL,
    description=DESCRIPTION,
    url=HOME_PAGE,
    project_urls={
        'Documentation': DOCUMENTATION,
        'Source': HOME_PAGE,
        'Tracker': f'{HOME_PAGE}/issues',
    },
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
        'thermodynamic-models', 'NRTL', 'UNIQUAC', 'SRK', 'PR', 'RK', 'vdW'
    ],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Education",
        "Programming Language :: Python :: 3.13",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ],
    python_requires='>=3.10',
)
