import jsmin
import minify_html
import oapackage
from oapackage.markup import oneliner as e


def formatArrayHyperlink(txt, lnk, filename):
    """Create HTML link to .oa file

    Args:
        txt (str): text to display
        lnk (str): url
        filename (str): .oa file on disk
    Returns:
        str: generated html

    """
    if oapackage.oahelper.oaIsBinary(filename):
        ss = e.a(txt + e.small(" (binary)"), href=lnk, class_="binarylink")
    else:
        ss = e.a(txt, href=lnk, class_="normal")
    return ss


def minifyJS(javascript_code: str) -> str:
    """Minify Javascript code"""
    minified = jsmin.jsmin(javascript_code)
    return minified


def minifyHTML(htmlcode: str, verbose: int = 0) -> str:
    """Minify html code"""
    htmlcode_minified = minify_html.minify(htmlcode, minify_js=True, minify_css=True)
    if verbose:
        print("minify: length %d -> %d" % (len(htmlcode), len(htmlcode_minified)))
    return htmlcode_minified
