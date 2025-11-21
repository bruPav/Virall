#!/usr/bin/env python3
"""
Virall - A comprehensive tool for viral genome analysis
including assembly, classification, gene prediction, and annotation.
"""

from setuptools import setup, find_packages
import re

# Read version from virall/__init__.py (single source of truth)
def get_version():
    with open("virall/__init__.py", "r", encoding="utf-8") as f:
        content = f.read()
        match = re.search(r'__version__\s*=\s*["\']([^"\']+)["\']', content)
        if match:
            return match.group(1)
        raise RuntimeError("Unable to find version string")

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="virall",
    version=get_version(),
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
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.11,<3.12",
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
