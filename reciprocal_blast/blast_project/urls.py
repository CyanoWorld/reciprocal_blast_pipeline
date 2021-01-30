from django.urls import path

from . import views

urlpatterns = [
    path('register',views.registrationPage,name='register'),
    path('login', views.loginPage, name='login'),
    path('logout', views.logoutUser, name='logout'),
    path('success',views.success_view,name='success_view'),
    path('<int:project_id>/pipeline_dashboard',views.pipeline_dashboard,name="pipeline_dashboard"),
    path('<int:project_id>/pipeline_nr_dashboard', views.pipeline_nr_dashboard, name="pipeline_nr_dashboard"),
    path('execute_snakefile',views.execute_snakefile,name='execute_snakefile'),
    path('execute_nr_snakefile', views.execute_nr_snakefile, name='execute_nr_snakefile'),
    path('project_creation',views.create_project,name='project_creation'),
    path('project_creation_based_on_nr_database', views.create_nr_based_project,name='project_creation_based_on_nr_database'),
    path('project_creation_based_on_uploaded_files', views.create_uploaded_based_project,name='project_creation_based_on_uploaded_files'),
    path('<int:project_id>/', views.project_details, name='project_details'),
    path('<int:project_id>/reciprocal_results', views.display_reciprocal_result_table, name='reciprocal_results'),
    path('<int:project_id>/delete',views.delete_project,name='delete_project'),
    path('failure',views.failure_view,name='failure_view'),
    path('species_taxid',views.species_taxid,name='species_taxid'),
    path('upload_databases',views.upload_databases,name='upload_databases'),
    path('', views.main, name='main'),
]