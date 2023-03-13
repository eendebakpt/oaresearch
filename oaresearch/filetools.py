import os
import shutil

import oapackage
from oapackage.oahelper import checkFiles


def copyOAfile(source, targetdir, target0, convert=None, zipfile=None, verbose=1, cache=1):
    """Copy an OA file, depending on arguments convert or compress the file

    Args:
        source (str): source file
        target (str): target output file
        convert (None, int or str): if None copy file. If a number, then convert to either text or binary. If a string, then 'T', 'B', or 'D'
        zipfile (None or bool): of True, compress the file with gzip
    Returns:
        target0final (str): output filename (without .gz extension)
    """
    if convert is None:
        targetfile = os.path.join(targetdir, target0)
        target0final = target0
        if not checkFiles(targetfile, cache):
            if verbose:
                print(f"copyfile {source} -> {targetfile}")
            shutil.copyfile(source, targetfile)
    else:
        na = oapackage.nArrayFile(source)
        if not (isinstance(convert, bytes) or isinstance(convert, str)):
            if na < convert:
                convert = "T"
            else:
                convert = "B"
        if zipfile is None:
            zipfile = convert == "B"
        else:
            if zipfile is not False:
                zipfile = na >= zipfile
        if not (
            convert == "TEXT"
            or convert == "BINARY"
            or convert == "B"
            or convert == "T"
            or convert == "D"
            or convert == "Z"
            or convert == "DIFF"
        ):
            raise NameError("copyOAfile: convert: should be T, B or D")
        if verbose >= 3:
            print(f"target0: {target0}, zipfile {zipfile}")
        if zipfile:
            if verbose:
                print("copyOAfile: converting to format %s" % convert)
            if target0.endswith(".gz"):
                target0final = target0
                target0 = target0final[:-3]
            else:
                target0final = target0 + ".gz"
            targetfilefinal = os.path.join(targetdir, target0final)
            targetfile = os.path.join(targetdir, target0)
        else:
            target0final = target0
            if target0final.endswith(".gz"):
                raise Exception("error: target file ends with .gz")
            targetfile = os.path.join(targetdir, target0)
            targetfilefinal = os.path.join(targetdir, target0)
        if verbose:
            print(f"copyOAfile: target0 {target0} -> {target0final} ")
        if verbose >= 2:
            print("copyOAfile: converting %s to %s (%d arrays, zip %d)" % (source, targetfile, na, zipfile))
        if checkFiles(targetfilefinal, cache):
            print("  copyOAfile: target file already exist")
        else:
            file_mode = oapackage.arrayfile_t_parseModeString(convert)
            oapackage.convert_array_file(source, targetfile, file_mode)
            if verbose >= 2:
                print("cmd: oaconvert -v 0 -f %s %s %s" % (convert, source, targetfile))
            if zipfile:
                os.system("gzip -f %s" % targetfile)
    return target0final
