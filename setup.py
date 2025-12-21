"""
PlastidTaxa: Scalable nucleid sequence taxonomy pipeline
Supports plastid (23S, 16S) and other nucleid targets
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="plastidtaxa",
    version="0.2.0",
    author="xaoscience",
    description="Nucleid-agnostic sequence taxonomy pipeline (plastid, mitochondrial, etc.)",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/xaoscience/PlastidTaxa",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.10",
    install_requires=[
        "numpy>=2.0",
        "scipy>=1.15",
        "pandas>=2.0",
        "matplotlib>=3.8",
        "seaborn>=0.13",
        "plotly>=6.0",
        "biopython>=1.87",  # XXE vulnerability fixed in 1.87+
        "defusedxml>=0.0.1",  # Mitigation for XXE in Bio.Entrez
        "scikit-learn>=1.0",
        "dada2>=1.30",  # Python port of DADA2
        "Click>=8.0",  # CLI framework
    ],
    extras_require={
        "dev": [
            "pytest>=7.0",
            "pytest-cov>=4.0",
            "black>=23.0",
            "flake8>=6.0",
            "sphinx>=5.0",
        ],
        "notebook": [
            "jupyter>=1.0",
            "notebook>=7.0",
            "ipython>=9.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "plastidtaxa=plastidtaxa.cli:main",
        ],
    },
)
