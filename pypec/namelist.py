"""
:mod:`pypec.namelist` -- IO for FORTRAN Namelists
=================================================

This module provides pythonic tools for reading and writing
fortran namelists.

Examples
--------

This module can be used to convert FORTRAN namelist files into
python ordered dictionary type objects. The individual namelists are
capitalized and their paramters are left as written. The actual
parameter is thus stored two levels deep in the dictionary.

>>> nl=read('examples/example_namelist.in')
>>> nl.keys()
['GPEC_INPUT', 'GPEC_CONTROL', 'GPEC_OUTPUT', 'GPEC_DIAGNOSE']
>>> nl['GPEC_CONTROL'].keys()
['resp_index', 'sing_spot', 'reg_flag', 'reg_spot', 'chebyshev_flag', 'nche']
>>> nl['GPEC_CONTROL']['reg_spot']
0.05

The parameters can be changed, added, or deleted as in any python
dictionary.

>>> nl['GPEC_CONTROL']['reg_spot'] = 0.01
>>> nl['GPEC_CONTROL']['reg_spot']
0.01
>>> del nl['GPEC_CONTROL']['nche']
>>> nl['GPEC_CONTROL'].keys()
['resp_index', 'sing_spot', 'reg_flag', 'reg_spot', 'chebyshev_flag']
>>> nl['GPEC_CONTROL']['nche'] = 30

The resulting namelists can then be re-written in ASCII format
for use with the FORTAN codes.

>>> write(nl,'examples/example_namelist_edit.in')
True

.. note::
 These examples can be tested by developers using ipython as follows:

   In [1]: import pypec.namelist,doctest

   In [2]: doctest.testmod(pypec.namelist,verbose=True)

"""

"""
    @package gpec
    @author NC Logan
    @email nikolas.logan@columbia.edu
"""

from collections import OrderedDict


class Objectify(object):
    """
    Class base to convert iterable to instance.

    """

    def __init__(self, d):
        """
        Recursively convert iterable to instance.

        :param d: iterable. Dictionary that will be converted.

        """
        for a, b in d.items():
            if isinstance(b, (list, tuple)):
                setattr(
                    self, a, [Objectify(x) if isinstance(x, dict) else x for x in b]
                )
            else:
                setattr(self, a, Objectify(b) if isinstance(b, dict) else b)

    def _todict(self):
        """
        Converts instance back to dictionary.

        """
        d = OrderedDict()
        for attr in dir(self):
            if not attr.startswith("__") and not attr == "_todict":
                val = getattr(self, attr)
                if isinstance(val, (list, tuple)):
                    d[attr] = [
                        x._todict() if isinstance(x, Objectify) else x for x in val
                    ]
                else:
                    d[attr] = val._todict() if isinstance(val, Objectify) else val
        return d


def _string_to_type(str):
    """
    Convert a string representing a fortran namelist value
    to the appropriate python type object.

    :param str  : str. String to convert

    :returns: obj. Python object of appropriate type.

    """
    try:
        pyobj = int(str)
    except ValueError:
        try:
            pyobj = float(str)
        except ValueError:
            if str.lower().startswith(".t") or str.lower().startswith("t"):
                pyobj = True
            elif str.lower().startswith(".f") or str.lower().startswith("f"):
                pyobj = False
            # complex written as tuple
            elif str.startswith("(") and str.endswith(")"):
                fpair = list(
                    map(
                        float,
                        str.translate(str.maketrans(dict.fromkeys("() "))).split(","),
                    )
                )
                pyobj = complex(*fpair)
            # lists written using spaces or repetition
            elif ("*" in str or " " in str or "," in str) and str.translate(
                str.maketrans(dict.fromkeys(",.* "))
            ).isalnum():
                pyobj = []
                str = str.replace(",", " ")
                for i in str.split():
                    if "*" in i:
                        n, v = i.split("*")
                        pyobj += int(n) * [_string_to_type(v)]
                    else:
                        pyobj += [_string_to_type(i)]
            # default to string
            else:
                pyobj = str.replace("'", "").replace('"', "")
    return pyobj


def read(filename, out_type="dict", comment="!", old=False):
    """
    Read in a file. Must have namelists in series, not embedded.
    Assumes name of list preceded by '&' and list ends with a single backslash.

    :param filename: str. Path to input namelist.
    :param out_type: str. Can be either 'dict' or 'obj'.
    :param comment: str. Lines starting with any of these characters are considered annotations.

    :returns: dict (object). Top keys (attributes) are namelists with sub-key(-attribute) parameter values.

    .. note:: The returned dict is actually "OrderedDict" type from
        collections module. If returning an object, the object is
        un-ordered.

    """

    if not out_type in ["dict", "obj"]:
        raise ValueError("Must output 'dict' or 'obj' type.")
    nldict = OrderedDict()
    if old:
        grp = "$"
        end = "/"  # $END gets mixed up with the $ group labeling
    else:
        grp = "&"
        end = "/"

    with open(filename, "r") as f:
        filestr = f.read()
    # handle $END block specification
    filestr = filestr.replace("$end", "$END").replace("$END", "/")
    # handle unnamed groups
    if grp not in filestr:
        print("No group names found: placing all in blank group.")
        filestr = grp + " \n" + filestr + "\n/"
    if "=" in filestr[0 : filestr.index(grp)]:
        print("Group with no name found in begining: placing in blank group.")
        filestr = (
            grp
            + " \n"
            + filestr[0 : filestr.index(grp)]
            + "\n"
            + end
            + "\n"
            + filestr[filestr.index(grp) :]
        )
    if "=" in filestr[filestr.rfind(end) :]:
        print("Group with no name found in end: placing in blank group.")
        filestr = (
            filestr[0 : filestr.rfind(end) + 1]
            + "\n"
            + grp
            + " \n"
            + filestr[filestr.rfind(end) + 1 :]
            + "\n"
            + end
        )
    # remove tabs and split lines
    nlist = (
        filestr.translate(str.maketrans(dict.fromkeys("\t")))
        .replace(grp, grp + "\n")
        .split("\n")
    )
    # prevent errors in searching for & or /, without condensing list values
    nlist = [nlistkey.strip() for nlistkey in nlist]

    start = 0
    stop = 0
    it = 0
    while stop < len(nlist):
        it += 1
        start = nlist.index(grp, stop) + 1
        stop = nlist.index(end, stop + 2) - 1
        name = nlist[start].strip().upper()
        nldict[name] = OrderedDict()
        # loop through each list, skiping empty lines
        # for item in [line for line in nlist[start+1:stop+1] if line.strip()]:
        for item in nlist[start + 1 : stop + 1]:
            if (item.lstrip() + " ")[0] in comment:  # comment
                iname, ival = item, None
                nldict[name][item] = None
            elif "=" in item:  # parameter
                iname = str.strip(item.split("=")[0])
                ival = item.split("=")[
                    1
                ]  # note there could be another '=' in the comments
                for c in comment:
                    ival = str.strip(ival.split(c)[0])  # allow same line comments
                try:
                    nldict[name][iname] = _string_to_type(ival)
                except:
                    raise ValueError("Failed to read " + name + " due to line: " + item)
            else:
                continue
        try:
            test = nlist.index(grp, stop)
        except ValueError:
            stop = len(nlist)
        if it > 100:
            stop = len(nlist)  # prevent inf loop if bug

    if out_type == "obj":
        return Objectify(nldict)
    return nldict


def write(nl, filename, old=False):
    """
    Write namelist object to file for fortran interface.

    :param nl: object. namelist object from read(). Can be dictionary or Objectify type class.
    :param filename: str. Path of file to be written.

    :returns: bool. True.

    """
    if type(nl) == Objectify:
        print("Warning: lists from objects are un-ordered.")
        nl = nl._todict()
    if old:
        start = "$"
        end = "$END"
    else:
        start = "&"
        end = "/"

    with open(filename, "w") as f:
        f.write("\n")  # vacuum code requires first line in vac.in
        for key in nl:
            f.write(start * bool(key) + key + "\n")
            for subkey in nl[key]:
                if type(nl[key][subkey]) == bool:
                    f.write(
                        "   " + subkey + " = ." + str(nl[key][subkey]).upper() + "."
                    )
                elif type(nl[key][subkey]) == complex:
                    fortranfmt = str(nl[key][subkey]).replace("+", ",").replace("j", "")
                    f.write("   " + subkey + " = " + fortranfmt + "")
                elif type(nl[key][subkey]) == str:
                    f.write("   " + subkey + " = '" + str(nl[key][subkey]) + "'")
                elif type(nl[key][subkey]) == list:
                    spacedlist = str(nl[key][subkey]).translate(
                        str.maketrans(dict.fromkeys(",[]"))
                    )
                    f.write("   " + subkey + " = " + spacedlist + "")
                elif type(nl[key][subkey]) is None:  # comments
                    f.write("   " + subkey)
                else:
                    f.write("   " + subkey + " = " + str(nl[key][subkey]))
                f.write("\n")
            f.write(end * bool(key) + "\n")

    return True


def _write_rst(nl, name="", filename="namelist.rst"):
    """
    Write namelist object to rst file for sphinx documentation
    including tooltips.

    :param nl       : object.  namelist object from read.
    :param name  : str. Name of namelist file.
    :param filename : str. Path of file to be written.

    :returns: bool. True.

    """
    import _tooltips_

    reload(_tooltips_)
    if type(nl) == Objectify:
        print("Warning: lists from objects are un-ordered.")
        nl = nl._todict()

    with open(filename, "w") as f:
        f.write(name.upper() + " Namelist Inputs\n")  # Title
        f.write("*" * len(name) + "****************\n\n")
        for key in nl:
            f.write(key + "\n")
            f.write("=" * len(key) + "\n\n")
            for subkey in nl[key]:
                if type(nl[key][subkey]) == bool:
                    f.write(
                        "**"
                        + subkey
                        + "**"
                        + " = ."
                        + str(nl[key][subkey]).upper()
                        + ".\n"
                    )
                elif type(nl[key][subkey]) == complex:
                    fortranfmt = str(nl[key][subkey]).replace("+", ",").replace("j", "")
                    f.write("**" + subkey + "**" + " = " + fortranfmt + "\n")
                elif type(nl[key][subkey]) == str:
                    f.write(
                        "**" + subkey + "**" + ' = "' + str(nl[key][subkey]) + '"\n'
                    )
                elif type(nl[key][subkey]) == list:
                    spacedlist = str(nl[key][subkey]).translate(
                        str.maketrans(dict.fromkeys(",[]"))
                    )
                    f.write("**" + subkey + "**" + " = " + spacedlist + "\n")
                else:
                    f.write("**" + subkey + "**" + " = " + str(nl[key][subkey]) + "\n")
                if subkey in _tooltips_.alltips:
                    f.write(
                        "  " + _tooltips_.alltips[subkey] + "\n"
                    )  # .replace('\n','')+'\n')
                else:
                    f.write("\n")
                f.write("\n")
            f.write("\n\n")

    return True
