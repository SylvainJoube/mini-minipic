"""Lib execution from command line."""

import sys

from libminipic.exceptions import MiniPICError
from libminipic.run import run
from libminipic.validate import validate


def execute(function):
    """Execute a function and catch usual errors."""

    try:
        function()
        return 0

    except MiniPICError as exception:
        print(exception)
        return 1

    except:
        print("Unexpected error")
        raise


def execute_validate():
    sys.exit(execute(validate))


def execute_run():
    sys.exit(execute(run))
