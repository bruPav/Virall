"""
Viral Genome Assembler - A comprehensive tool for viral genome assembly
from short and long sequencing reads with high accuracy.
"""

__version__ = "0.2.2"
__author__ = "Viral Assembly Team"
__email__ = "bruno.pavletic@gmail.com"

from .core.assembler import ViralAssembler
from .core.preprocessor import Preprocessor
from .core.viral_identifier import ViralIdentifier
from .core.validator import AssemblyValidator

__all__ = [
    "ViralAssembler",
    "Preprocessor", 
    "ViralIdentifier",
    "AssemblyValidator"
]
