from django.urls import path
# Import views from the current application package.
from . import views

urlpatterns = [
    # The site's home page handles the form and all three analysis actions.
    # The name can be used with Django's ``url`` template tag or reverse().
    path('', views.predict_go, name='predict_go'),
]
