"""Load unmodified upstream components without its inference-only imports."""

import importlib
import sys
import types
from pathlib import Path


def load_reference():
    package_name = "_polii_polii_reference"
    location = Path(__file__).resolve().parents[1] / "upstream" / "polii"
    if not (location / "model.py").is_file():
        raise FileNotFoundError(f"Missing reference PolII source: {location}")
    if package_name not in sys.modules:
        package = types.ModuleType(package_name)
        package.__path__ = [str(location)]
        sys.modules[package_name] = package
    return importlib.import_module(f"{package_name}.model")
