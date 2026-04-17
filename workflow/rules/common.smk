# workflow/rules/common.smk
# Shared helpers imported by the main Snakefile before any other rules.

import os
from pathlib import Path


def p(*parts) -> str:
    """Join path parts — avoids string literals starting with '/' in rule directives."""
    return str(Path(*parts))
