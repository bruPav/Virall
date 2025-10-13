"""
Viral Genome Assembler - A comprehensive tool for viral genome assembly
from short and long sequencing reads with high accuracy.
"""

__version__ = "1.0.0"
__author__ = "Viral Assembly Team"
__email__ = "team@viralassembler.com"

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
