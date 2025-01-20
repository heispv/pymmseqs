# pymseqs/__init__.py

from .binary import get_mmseqs_binary
from .runner import run_mmseqs_command

__all__ = [
    "get_mmseqs_binary",
    "run_mmseqs_command",
]
