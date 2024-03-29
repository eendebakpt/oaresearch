""" Misc functions

Pieter Eendebak <pieter.eendebak@gmail.com>

"""


import numpy as np

try:
    pass
except BaseException:
    pass

import oapackage
import oapackage.markup as markup

# %%


def array2cpp(array: oapackage.array_link, padding: str = "   "):
    """Convert array to C++ initialization code"""
    ss = padding + "array_link array (%d, %d, 0);" % (array.n_rows, array.n_columns) + " \n"
    ss += padding + "int array_data_tmp[] = {%s};\n" % (",".join(["%d" % v for v in np.array(array).T.flatten()]))
    ss += padding + "array.setarraydata (array_data_tmp, array.n_rows * array.n_columns);\n"
    return ss


def array2html(
    X, header=1, tablestyle="border-collapse: collapse;", trclass="", tdstyle="", trstyle="", thstyle="", comment=None
):
    """Convert Numpy array to HTML table

    Arguments
    ---------
        X : numpy array
            array to be converted
        header : integer
            use header or not
    Returns
    -------
        page : markup html object
            generated table in HTML

    """
    page = markup.page()
    page.add("<!-- Created by array2html -->\n")
    if comment is not None:
        page.add("<!-- %s-->\n" % comment)

    page.table(style=tablestyle)
    offset = 0
    nc = X.shape[1]
    nr = X.shape[0]

    if isinstance(trstyle, str):
        trstyle = [trstyle] * nr
    if isinstance(trclass, str):
        trclass = [trclass] * nr

    ri = 0
    if header:
        page.tr(style="font-weight: bold; border-bottom: solid 1px black;" + trstyle[ri], class_=trclass[ri])
        ri = ri + 1
        for ii in range(nc):
            if isinstance(X[offset, ii], tuple):
                print("array2html: tuple instance")
                page.th(X[offset, ii][0], style=thstyle + X[offset, ii][1])
            else:
                page.th(X[offset, ii], style=thstyle)
        page.tr.close()
        offset = offset + 1

    nr = X.shape[0] - offset
    for r in range(nr):
        page.tr(style=trstyle[ri], _class=trclass[ri])
        for ii in range(nc):
            if isinstance(X[offset, ii], tuple):
                page.td(X[offset, ii][0], style=tdstyle + X[offset, ii][1])
            else:
                page.td(X[offset, ii], style=tdstyle)

        page.tr.close()
        offset = offset + 1
        ri = ri + 1
    page.table.close()
    return page


# %%


def citation(paper, style="brief"):
    """Return citation in html format
    Args:
        paper (str): paper to be cited
        style (str): brief or full (with authors and journal)
    """
    if paper == "complete":
        return markup.oneliner.a(
            "Complete Enumeration of Pure-Level and Mixed-Level Orthogonal Arrays",
            href="https://doi.org/10.1002/jcd.20236",
        )
    elif paper == "conference" or paper == "cisomorphism":
        if style == "full":
            hyperlink = markup.oneliner.a(
                "A Classification Criterion for Definitive Screening Designs", href="https://doi.org/10.1214/18-AOS1723"
            )
            return hyperlink + ", E.D. Schoen, P.T. Eendebak, P. Goos, Annals of Statistics, 2019"
        else:
            return markup.oneliner.a(
                "A Classification Criterion for Definitive Screening Designs", href="https://doi.org/10.1214/18-AOS1723"
            )
    elif paper == "conference enumeration" or paper == "cenumeration":
        if style == "full":
            return (
                markup.oneliner.a("Enumeration and Classification of Definitive Screening Designs", href="")
                + ", Eric D. Schoen, Pieter T. Eendebak, Alan Vazquez-Alcocer, Peter Goos, in preparation"
            )
        else:
            return "Enumeration and Classification of Definitive Screening Designs, in preparation"
            return markup.oneliner.a("Enumeration and Classification of Definitive Screening Designs", href="")
    else:
        raise Exception("paper not known")


def oaCssStyle(addframe=False):
    ss = """
    /* Style file for Orthognal Arrays page
 * Pieter Eendebak <pieter.eendebak@gmail.com>
 *
 */

body {font-family: Helvetica,Arial,sans-serif; }
a { color: #0000cc; text-decoration: none;}
a:hover { color: #110149; text-decoration: none;}
span.series a { color: blue; text-decoration: none;}
a.binarylink { color: #7722FF; }


.top th {padding-bottom: .4em; }
tr.odd {}
//tr.block{ background-color: red; }
tr.block { border-top: 1px #000 solid; /* top border only */
border-color: #e0e0e0; }
.block td { padding-top: .4em;}
//tr.even {background-color: #fefefe; }
tr.even {background-color: #fafafa; }

"""
    if addframe:
        ss = '<style type="text/css">\n' + ss + "</style>\n"
        return ss
    else:
        return ss
