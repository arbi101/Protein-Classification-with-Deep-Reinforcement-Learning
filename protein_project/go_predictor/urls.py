from django.urls import path
from . import views

urlpatterns = [
    path('', views.predict_go, name='predict_go'),
    path('structure/', views.structure_2d, name='structure_2d'),
]