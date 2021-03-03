'''
Content: URLs and associated view functions.
For further information see Django developer guidelines
'''

from django.urls import path
from . import views

urlpatterns = [
    path('register',views.registration_view,name='register'),
    path('login', views.login_view, name='login'),
    path('logout', views.logout_user, name='logout'),
    path('success',views.success_view,name='success_view'),
    path('<int:project_id>/pipeline_dashboard',views.pipeline_dashboard,name="pipeline_dashboard"),
    path('<int:project_id>/pipeline_nr_dashboard', views.pipeline_nr_dashboard, name="pipeline_nr_dashboard"),
    path('execute_snakefile',views.execute_snakefile,name='execute_snakefile'),
    path('execute_nr_snakefile', views.execute_nr_snakefile, name='execute_nr_snakefile'),
    path('<int:project_id>/', views.project_details, name='project_details'),
    path('<int:project_id>/reciprocal_results', views.display_reciprocal_result_table, name='reciprocal_results'),
    path('<int:project_id>/delete',views.delete_project,name='delete_project'),
    path('failure',views.failure_view,name='failure_view'),
    path('species_taxid',views.species_taxid,name='species_taxid'),
    path('upload_databases',views.upload_databases,name='upload_databases'),
    path('project_creation',views.project_creation,name='project_creation'),
    path('refseq_database_transactions',views.refseq_database_download,name='refseq_database_transactions'),
    path('download_refseq_summary',views.download_refseq_file,name='download_refseq_summary'),
    path('delete_refseq_summary_file',views.delete_refseq_summary_file,name='delete_refseq_summary_file'),
    path('', views.main, name='main'),
]