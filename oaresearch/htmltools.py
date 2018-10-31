import oapackage
from oapackage.markup import oneliner as e


def formatArrayHyperlink(txt, lnk, filename):
    """ Create HTML link to .oa file

    Args:
        txt (str): text to display
        lnk (str): url
        filename (str): .oa file on disk
    Returns:
        str: generated html

    """
    if oapackage.oahelper.oaIsBinary(filename):
        ss = e.a(txt + e.small(' (binary)'), href=lnk, class_='binarylink')
    else:
        ss = e.a(txt, href=lnk, class_='normal')
    return ss
