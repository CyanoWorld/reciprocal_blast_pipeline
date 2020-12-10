from django.urls import path

from . import views

urlpatterns = [
    path('register',views.registrationPage,name='register'),
    path('login', views.loginPage, name='login'),
    path('logout', views.logoutUser, name='logout'),
    path('success',views.success_view,name='success_view'),
    path('<int:project_id>/pipeline_dashboard',views.pipeline_dashboard,name="pipeline_dashboard"),
    path('build_snakefile', views.build_snakefile,name='build_snakefile'),
    path('remove_snakefile',views.remove_snakefile,name='remove_snakefile'),
    path('execute_snakefile',views.execute_snakefile,name='execute_snakefile'),
    path('project_creation',views.create_project,name='project_creation'),
    path('<int:project_id>/', views.project_details, name='project_details'),
    path('<int:project_id>/delete',views.delete_project,name='delete_project'),
    path('failure',views.failure_view,name='failure_view'),
    path('', views.main, name='main'),
]