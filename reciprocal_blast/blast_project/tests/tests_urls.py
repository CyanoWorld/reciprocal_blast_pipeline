from django.test import SimpleTestCase
from django.urls import resolve, reverse
#no error
from blast_project.views import main, project_details, project_creation


class TestUrls(SimpleTestCase):

    def test_main_url_resolves(self):
        url = reverse('main')
        print(resolve(url))
        self.assertEquals(resolve(url).func, main)

    def test_project_details_url_resolves(self):
        url = reverse('project_details',args='1')
        print(resolve(url))
        self.assertEquals(resolve(url).func,project_details)

    def test_create_project_url_resolves(self):
        url = reverse('project_creation')
        print(resolve(url))
        self.assertEquals(resolve(url).func,project_creation)