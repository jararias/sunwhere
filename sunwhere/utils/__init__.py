
import importlib


def safe_import(name, package=None):
    try:
        module = importlib.import_module(name, package)
    except (ImportError, ModuleNotFoundError):
        module = None
    return module
