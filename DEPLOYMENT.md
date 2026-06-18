# Deployment Notes

This is a Django application, so deploy it as a Python web service.

## Render

Create a new Web Service from this repository and use these settings:

```text
Build Command:
pip install -r requirements.txt && cd protein_project && python manage.py collectstatic --noinput && python manage.py migrate

Start Command:
cd protein_project && gunicorn protein_project.wsgi:application --workers 1 --timeout 120
```

Add this environment variable in Render:

```text
SECRET_KEY=<generate-a-long-secret-key>
```

Optional environment variables:

```text
DEBUG=False
ALLOWED_HOSTS=your-custom-domain.com
CSRF_TRUSTED_ORIGINS=https://your-custom-domain.com
TORCH_NUM_THREADS=1
DQN_ROLLOUTS=20
MAX_DQN_ROLLOUTS=40
MAX_STRUCTURE_SEQUENCES=1
```

For a quick secret key:

```bash
cd protein_project
python manage.py shell -c "from django.core.management.utils import get_random_secret_key; print(get_random_secret_key())"
```

The app currently uses SQLite. That is fine for a demo, but persistent production data should use a managed database.
