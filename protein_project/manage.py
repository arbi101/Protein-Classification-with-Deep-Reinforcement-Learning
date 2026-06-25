#!/usr/bin/env python
"""Django's command-line utility for administrative tasks."""
import os
import sys


def main():
    """Run administrative tasks."""
    # Tell Django which settings module must be loaded. ``setdefault`` keeps
    # any value explicitly provided by the environment (for example in tests).
    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'protein_project.settings')
    try:
        # This function parses commands such as runserver, migrate and test.
        from django.core.management import execute_from_command_line
    except ImportError as exc:
        # Replace a cryptic import failure with guidance about Django and the
        # virtual environment, while preserving the original exception.
        raise ImportError(
            "Couldn't import Django. Are you sure it's installed and "
            "available on your PYTHONPATH environment variable? Did you "
            "forget to activate a virtual environment?"
        ) from exc
    # Forward the command-line arguments to Django's command dispatcher.
    execute_from_command_line(sys.argv)


if __name__ == '__main__':
    # Run only when this file is executed directly, not when it is imported.
    main()
