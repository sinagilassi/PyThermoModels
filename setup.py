from setuptools import setup, find_packages
import codecs
import os

here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as fh:
    long_description = "\n" + fh.read()

APP_NAME = 'PyThermoModels'
VERSION = '1.0.0a1'
DESCRIPTION = 'PyThermoModels is a lightweight and user-friendly Python package for thermodynamic modeling.'
LONG_DESCRIPTION = "PyThermoModels is an open-source Python package designed to facilitate thermodynamic modeling and calculations. This package provides a comprehensive and user-friendly interface to popular thermodynamic models, enabling quick and accurate estimation of key properties."

# Setting up
setup(
    name=APP_NAME,
    version=VERSION,
    author="Sina Gilassi",
    author_email="<sina.gilassi@gmail.com>",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=long_description,
    packages=find_packages(exclude=['tests', '*.tests', '*.tests.*']),
    include_package_data=True,  # Make sure to include non-Python files
    # Add both config and data files
    package_data={'': ['config/*.yml', 'config/*.csv', 'data/*.csv',
                       'data/*.yml', 'plugin/*.yml', 'plugin/*.csv']},
    license='MIT',
    install_requires=['pandas', 'pillow', 'requests',
                      'urllib3', 'matplotlib', 'numpy', 'PyYAML', 'sympy', 'PyThermoDB'],
    keywords=['python', 'chemical engineering', 'thermodynamics',
              'PyThermoDB'],
    classifiers=[
        "Development Status :: 1 - Alpha",
        "Intended Audience :: Education",
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ],
    python_requires='>=3.6',
)
