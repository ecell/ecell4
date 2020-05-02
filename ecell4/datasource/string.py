import urllib
import time

import requests


class StringDataSource:
    """
    See http://string-db.org/cgi/help.pl

    """

    wait = 1.0

    def __init__(self):
        pass

    @classmethod
    def network(cls, entity, output_format='highres_image', **optional_parameters):
        """
        Display STRING network image.

        Parameters
        ----------
        output_format : str, optional
            'image' or 'highres_image'.

        """
        if isinstance(entity, str):
            identifiers = [entity]
        else:
            raise TypeError(f'The given identifier [{entity}] is invalid.')

        parameters = {'identifiers': '\r'.join(identifiers)}
        parameters.update(optional_parameters)

        url = "https://string-db.org/api/{}/network?{}".format(
                output_format, urllib.parse.urlencode(parameters))
        response = requests.get(url)
        time.sleep(cls.wait)

        from IPython.display import Image, display_png
        display_png(Image(response.content))

    @classmethod
    def interaction_partners(cls, entity, **optional_parameters):
        if isinstance(entity, str):
            identifiers = [entity]
        else:
            raise TypeError(f'The given identifier [{entity}] is invalid.')

        parameters = {'identifiers': '\r'.join(identifiers)}
        parameters.update(optional_parameters)

        output_format = 'json'
        url = "https://string-db.org/api/{}/interaction_partners?{}".format(
                output_format, urllib.parse.urlencode(parameters))
        response = requests.get(url)
        time.sleep(cls.wait)
        return response.json()
