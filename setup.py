#!/usr/bin/env python3
"""
Virall - A comprehensive tool for viral genome analysis
including assembly, classification, gene prediction, and annotation.
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="virall",
    version="0.1.0",
    author="Virall Team",
    author_email="bruno.pavletic@gmail.com, britaniadiazf@gmail.com",
    description="A comprehensive tool for viral genome analysis including assembly, classification, gene prediction, and annotation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bruPav/virall",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "virall=virall.cli:main",
        ],
    },
    include_package_data=True,
    package_data={
        "virall": ["data/*", "configs/*"],
    },
)
