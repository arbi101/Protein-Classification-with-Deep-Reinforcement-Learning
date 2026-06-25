"""
ASGI config for protein_project project.

It exposes the ASGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/6.0/howto/deployment/asgi/
"""

import os

from django.core.asgi import get_asgi_application

# Select this project's configuration before creating the server callable.
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'protein_project.settings')

# ASGI entry point for asynchronous-capable servers such as Uvicorn.
application = get_asgi_application()
