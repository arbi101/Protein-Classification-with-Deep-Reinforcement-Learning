"""
WSGI config for protein_project project.

It exposes the WSGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/6.0/howto/deployment/wsgi/
"""

import os

from django.core.wsgi import get_wsgi_application

# Select this project's configuration before creating the server callable.
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'protein_project.settings')

# WSGI entry point used by the Gunicorn production server.
application = get_wsgi_application()
