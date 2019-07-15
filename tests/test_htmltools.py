# -*- coding: utf-8 -*-
"""
Pieter Eendebak <pieter.eendebak@gmail.com>

"""

# %%
import unittest
from unittest.mock import patch
import io

import oaresearch.htmltools

# %%


class TestHTMLtools(unittest.TestCase):

    def test_formatArrayHyperlink(self):
        txt = 'hello'
        url = 'http://www.world'
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            html = oaresearch.htmltools.formatArrayHyperlink(
                txt, url, 'arrayfile')
            self.assertTrue(html.startswith('<a '))
            self.assertTrue(txt in html)
            self.assertTrue(url in html)


if __name__ == '__main__':
    unittest.main()
