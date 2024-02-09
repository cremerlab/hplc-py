from setuptools import setup, find_packages
import pathlib

__version__ = "0.2.5"

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

setup(
    name="hplc-py",
    version=__version__,
    long_description=README,
    description="Python utilities for the processing and quantification of chromatograms from High Performance Liquid Chromatography (HPLC).",
    long_description_content_type='text/markdown',
    url="https://github.com/cremerlab/hplc-py",
    license="GPLv3",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
    ],
    author="Griffin Chure",
    author_email="griffinchure@gmail.com",
    packages=find_packages(
        exclude=('docs', 'doc', 'sandbox', 'dev', 'hplc.egg-info')),
    include_package_data=True,
    install_requires=[
        "matplotlib>=3.7.0",
        "numpy>=1.24.4",
        "pandas>=2.0.3",
        "scipy>=1.10.0",
        "seaborn>=0.12.2",
        "termcolor>=2.3.0",
        "tqdm>=4.64.1",
    ],
)
