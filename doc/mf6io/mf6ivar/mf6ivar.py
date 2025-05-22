# list of input issues
#
# general
#   auxiliary should be AUXNAMES auxnames(naux)
#   differentiate somehow a block that has repeating records from other block types?
#   should mark time series variables in some way special.  Maybe an *?
#   should we begin all character string variables with a 'c'?
#
# mfsim.nam
#   add dimensions block to sim name file? (NMODELS, NEXCHANGES, NSOLUTIONGROUPS)
#   add NSLNMODELS TO NUMERICAL line so we now how many models to read
#   mxiter is inside solution_group (label the block instead?)
#   change 'NUMERICAL' to 'SMS'
#
# tdis
#   read perlen, nstp, tsmult separately as arrays
#
# gwfgwf
#   change blockname EXCHANGEDATA? DATA or LIST or GWFGWFDATA
#   change auxiliary to auxnames(naux)
#   Move to different part of user guide?
#
# sms
#
# dis
#
# disu
#
# ic
#   no options, so what do we do?  allow empty options blocks?  dummy option?
#
# chd
#   change to OBS8_FILENAME obs8_filename
#   change AUXILIARY to AUXNAMES
#
# wel
#   change to OBS8_FILENAME obs8_filename
#   change AUXILIARY to AUXNAMES
#
# drn
#   change to OBS8_FILENAME obs8_filename
#   change AUXILIARY to AUXNAMES
#
# riv
#   change to OBS8_FILENAME obs8_filename
#   change AUXILIARY to AUXNAMES
#
# ghb
#   change to OBS8_FILENAME obs8_filename
#   change AUXILIARY to AUXNAMES
#
# rch
#   change to OBS8_FILENAME obs8_filename
#   change AUXILIARY to AUXNAMES
#   do not support array-based input for DISU
#   looks like code can be cleaned up quite a bit if we require lists for DISU
#   for array-based input, should arrays be optional?
#
# evt
#   change to OBS8_FILENAME obs8_filename
#   change AUXILIARY to AUXNAMES
#   do not support array-based input for DISU
#   looks like code can be cleaned up quite a bit if we require lists for DISU
#   for array-based input, should arrays be optional?
#   do not support segmented ET for array-based input
#
# maw
#   need to remove the word WELL as the first thing in the WELLS block
#   change to OBS8_FILENAME obs8_filename
#   change AUXILIARY to AUXNAMES
#   suggest renaming ngwfnodes to ncon or nconn
#   change 'WELLS' to 'MAWDATA' (do this consistently throughout SFR, LAK, etc.?)
#   change 'WELL_CONNECTIONS' to 'CONNECTIONS'
#   eliminate 'STEADY-STATE' keyword from period block
#   Change ACTIVE, INACTIVE, and CONSTANT into values for STATUS.
#
# sfr
#   change to OBS8_FILENAME obs8_filename
#   change AUXILIARY to AUXNAMES
#   change from SFR_OUTPUT to SFROUTPUT_FILENAME?  or OUTPUT_FILENAME
#   rno does not correspond to implicit integer definition
#   unit conversion used in sfr, length/time conversion used in lake
#   need to implement STATUS
#
# lak
#   change to OBS8_FILENAME obs8_filename
#   change AUXILIARY to AUXNAMES
#   unit conversion used in sfr, length/time conversion used in lake
#   change time_conversion to timefactor?
#   change LAKE_CONNECTIONS to CONNECTIONS
#   change LAKE_TABLES to TABLES
#   change FILE ctabname to TABLE_FILENAME table_filename
#   STATUS not implemented yet, but it is described in input instructions
#   invert indicates an integer variable.  change to dinvert?
#   time series variables are listed as "real or character ..."
#       should just be double precision
#   capitalize example input file words that are recognized by mf6
#
# uzf
#   change to OBS8_FILENAME obs8_filename
#   change AUXILIARY to AUXNAMES
#   Remove steady-state / transient flag
#   Change DATA block to UZFDATA
#   boundname not implemented.  Should be implemented for uzfdata?
#   aux not implemented.  should be implemented for period block?
#   implement uzfsettings approach?  Or stick with a simple list?
#   example shows 'uzf' keyword as first item in period block
#   combine SIMULATE_ET, LINEAR_GWET, and SQUARE_GWET into a single line?
#
# mvr
#   change maxpackages to npackages
#   Included WEL, DRN, RIV, GHB as providers,
#       though that is not supported in the code yet
#
# oc
#   output control rewritten entirely, and implemented in the code
#

# DEFINITION FILE KEYWORDS
# block :: name of block
# name :: variable name
# in_record :: optional True or False, False if not specified
# type :: recarray, record, keyword, integer, double precision, keystring
# tagged :: optional True or False, True if not specified.
#           If tagged, then keyword comes before value
# shape :: (size), optional, only required for arrays
# valid :: description of valid values
# reader :: urword, readarray, u1dint, ...
# optional :: optional True or False, False if not specified
# longname :: long name for variable
# description :: description for variable, REPLACE tag indicates that
#                description will come from common.dfn


import os
import re
import shutil
import textwrap
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from collections import OrderedDict
from pathlib import Path


def parse_mf6var_file(fname):
    lines = open(fname, "r").readlines()
    vardict = OrderedDict()
    vd = {}

    for line in lines:
        # skip blank lines
        if len(line.strip()) == 0:
            if len(vd) > 0:
                name = vd["name"]
                if "block" in vd:
                    block = vd["block"]
                    key = (name, block)
                else:
                    key = name
                if key in vardict:
                    raise ValueError(f"Variable already exists in dictionary: {key}")
                vardict[key] = vd
            vd = {}
            continue

        # skip comments
        if "#" in line.strip()[0]:
            continue

        ll = line.strip().split()
        if len(ll) > 1:
            k = ll[0]
            istart = line.index(" ")
            v = line[istart:].strip()
            if k in vd:
                raise ValueError(f"Attribute already exists in dictionary: {k}")
            vd[k] = v

    if len(vd) > 0:
        name = vd["name"]
        if "block" in vd:
            block = vd["block"]
            key = (name, block)
        else:
            key = name
        if key in vardict:
            raise ValueError(f"Variable already exists in dictionary: {k}")
        vardict[key] = vd
    return vardict


MF6IVAR_DIR_PATH = Path(__file__).parent
MF6IO_DIR_PATH = Path(__file__).parents[1]
DFNS_DIR_PATH = MF6IVAR_DIR_PATH / "dfn"
EXAMPLES_DIR_PATH = MF6IVAR_DIR_PATH / "examples"
MD_DIR_PATH = MF6IVAR_DIR_PATH / "md"
TEX_DIR_PATH = MF6IVAR_DIR_PATH / "tex"
RTD_DOC_DIR_PATH = Path(__file__).parents[3] / ".build_rtd_docs" / "_mf6io"
COMMON_DFN_PATH = parse_mf6var_file(DFNS_DIR_PATH / "common.dfn")
COMMON_DIR_PATH = MF6IVAR_DIR_PATH.parent.parent / "Common"
DEFAULT_MODELS = ["gwf", "gwt", "gwe", "prt"]
DEVELOP_MODELS = ["chf", "olf"]
if (MF6IO_DIR_PATH / "develop.version").is_file():
    DEFAULT_MODELS += DEVELOP_MODELS
VALID_TYPES = [
    "integer",
    "double precision",
    "string",
    "keystring",
    "keyword",
    "recarray",
    "record",
]

MD_DIR_PATH.mkdir(exist_ok=True)
TEX_DIR_PATH.mkdir(exist_ok=True)


def block_entry(varname, block, vardict, prefix="  "):
    key = (varname, block)
    v = vardict[key]

    s = f"{varname.upper()}"
    if "tagged" in v:
        if v["tagged"] == "false":
            s = ""

    # set up the time series marker @
    tsmarker = ""
    if "time_series" in v:
        if v["time_series"] == "true":
            tsmarker = "@"
    extmarker = ""
    if "extended" in v:
        if v["extended"] == "true":
            extmarker = "$"
    if "netcdf" in v:
        if v["netcdf"] == "true":
            extmarker = "$"

    # check valid type
    vtype = v["type"]
    if vtype == "double precision":
        pass
    elif " " in vtype:
        vtype = vtype.split(" ", 1)[0]
    if vtype not in VALID_TYPES:
        raise ValueError(f"{key}: {vtype!r} is not a valid type from {VALID_TYPES}")

    # record or recarray
    if v["type"].startswith("rec"):
        varnames = v["type"].strip().split()[1:]
        s = ""
        for vn in varnames:
            blockentry = block_entry(vn, block, vardict, prefix="")
            s += f"{blockentry.strip()} "
        if v["type"].startswith("recarray"):
            s = s.strip()
            s = f"{s}\n{prefix}{s}\n{prefix}..."

    # layered and netcdf
    elif v["reader"] in ["readarray", "u1ddbl", "u2ddbl", "u1dint"]:
        shape = v["shape"]
        reader = v["reader"].upper()
        layered = ""
        if "layered" in v:
            if v["layered"] == "true":
                layered = " [LAYERED]"
        if "netcdf" in v:
            if v["netcdf"] == "true":
                layered = layered + f" {extmarker}[NETCDF]{extmarker}"
        s = f"{s}{layered}\n{prefix}{prefix}<{varname}{shape}> -- {reader}"

    # timeseries, extended color annotation
    else:
        vtmp = varname
        if tsmarker != "" and v["type"] != "keyword":
            if "shape" in v:
                shape = v["shape"]
                vtmp += shape
            s = f"{s} <{tsmarker}{vtmp}{tsmarker}>"
        elif extmarker != "":
            if v["type"] != "keyword":
                if "shape" in v:
                    shape = v["shape"]
                    vtmp += shape
                s = f"{extmarker}{s}{extmarker} <{extmarker}{vtmp}{extmarker}>"
            else:
                s = f"{extmarker}{s}{extmarker}"
        elif v["type"] != "keyword":
            s = f"{s} <{vtmp}>"

    # if optional, wrap string in square brackets
    if "optional" in v:
        if v["optional"] == "true":
            s = f"[{s.strip()}]"

    # prepend with prefix and return string
    s = f"{prefix}{s}"
    return s


def write_block(vardict, block, blk_var_list, varexcludeprefix=None, indent=None):
    prepend = "" if indent is None else indent * " "
    s = prepend + f"BEGIN {block.upper()}"
    for variable in blk_var_list:
        ts = block_entry(variable[0], block, vardict).strip()
        if variable[1]:
            s = f"{s} [{ts}]"
        else:
            s = f"{s} {ts}"
    s += "\n"
    for key in vardict.keys():
        name, b = key
        v = vardict[key]
        if b == block:
            addv = True
            if varexcludeprefix is not None:
                # exclude variables that start with `dev_`.  These are
                # develop options that shouldn't go into the documentation.
                n = name.upper()
                if n.startswith(varexcludeprefix.upper()):
                    addv = False
            if v.get("in_record", "") == "true":
                # do not separately include this variable
                # because it is part of a record
                addv = False
            if v.get("block_variable", "") == "true":
                # do not separately include this variable
                # because it is part of a record
                addv = False
            if v.get("deprecated", "") != "" or v.get("removed", "") != "":
                addv = False
            if addv:
                ts = block_entry(name, block, vardict, prefix="  " + prepend)
                s += f"{ts}\n"
    s += prepend + f"END {block.upper()}"
    return s


def get_description(desc):
    """
    Check to see if the description is in common.dfn, and make the text
    substitutions if so.

    """
    if desc.strip().split()[0] == "REPLACE":
        bcoption = desc.strip().split()[1]
        constantstring = COMMON_DFN_PATH[bcoption]["description"]
        istart = desc.index("{")
        istop = desc.rfind("}") + 1
        d = eval(desc[istart:istop])
        for k in d:
            v = d[k]
            constantstring = constantstring.replace(k, v)
        desc = constantstring + desc[istop:]
    return desc


def write_desc(vardict, block, blk_var_list, varexcludeprefix=None):
    s = ""
    for name, b in vardict.keys():
        v = vardict[(name, b)]
        if v["block"] == block:
            if v.get("block_variable"):
                optional = "optional" in v and v["optional"] == "true"
                blk_var_list.append((v["name"], optional))
            addv = True
            if varexcludeprefix is not None:
                # exclude variables that start with `dev_`.  These are
                # develop options that shouldn't go into the documentation.
                n = name.upper()
                if n.startswith(varexcludeprefix.upper()):
                    addv = False
            if v["type"].startswith("rec"):
                addv = False
            if v.get("deprecated", "") != "":
                addv = False
            if v.get("removed", "") != "":
                addv = False
            if addv:
                if v["type"] == "keyword":
                    n = name.upper()
                else:
                    if "tagged" in v:
                        # could be used in future to write tag and name
                        n = f"{name}"
                    else:
                        n = name
                n = n.replace("_", "\\_")
                if "description" in v:
                    desc = get_description(v["description"])
                else:
                    msg = ""
                    for k, v in v.items():
                        msg += f"  {k}: {v}\n"
                    print(msg)
                    raise Exception(msg)
                ss = "\\texttt{" + n + "}---" + desc
                if "time_series" in v:
                    if v["time_series"] == "true":
                        fmt = "\\textcolor{blue}\{\}"
                        ss = "\\textcolor{blue}{" + ss + "}"
                        # \textcolor{declared-color}{text}
                if "extended" in v:
                    if v["extended"] == "true":
                        fmt = "\\textcolor{red}\{\}"
                        ss = "\\textcolor{red}{" + ss + "}"
                s += "\\item " + ss + "\n\n"

                t = v["type"]
                if t.startswith("keystring"):
                    # s += '\\begin{verbatim}\n'
                    s += "\\begin{lstlisting}[style=blockdefinition]\n"
                    for vn in t.strip().split()[1:]:
                        if (
                            "removed" in vardict[(vn, block)]
                            or "deprecated" in vardict[(vn, block)]
                        ):
                            continue
                        blockentry = block_entry(vn, block, vardict, "")
                        s += f"{blockentry}\n"
                    # s += '\\end{verbatim}\n\n'
                    s += "\\end{lstlisting}\n\n"

    return s


def write_desc_md(vardict, block, blk_var_list, varexcludeprefix=None):
    s = ""
    for name, b in vardict.keys():
        v = vardict[(name, b)]
        if v["block"] == block:
            if v.get("block_variable"):
                optional = "optional" in v and v["optional"] == "true"
                blk_var_list.append((v["name"], optional))
            addv = True
            if varexcludeprefix is not None:
                # exclude variables that start with `dev_`.  These are
                # develop options that shouldn't go into the documentation.
                n = name.upper()
                if n.startswith(varexcludeprefix.upper()):
                    addv = False
            if v["type"].startswith("rec"):
                addv = False
            if v.get("deprecated", "") != "":
                addv = False
            if v.get("removed", "") != "":
                addv = False
            if addv:
                if v["type"] == "keyword":
                    n = name.upper()
                else:
                    if "tagged" in v:
                        # could be used in future to write tag and name
                        n = f"{name}"
                    else:
                        n = name
                if "description" in v:
                    desc = get_description(v["description"])
                else:
                    msg = ""
                    for k, v in v.items():
                        msg += f"  {k}: {v}\n"
                    print(msg)
                    raise Exception(msg)
                desc = md_replace(desc)
                ss = "`" + n + "` " + desc
                if "time_series" in v:
                    if v["time_series"] == "true":
                        ss = '<span style="color:blue">' + ss + "</span>"
                if "extended" in v:
                    if v["extended"] == "true":
                        ss = '<span style="color:red">' + ss + "</span>"
                s += "  * " + ss + "\n\n"

                t = v["type"]
                if t.startswith("keystring"):
                    for vn in t.strip().split()[1:]:
                        blockentry = md_replace(
                            block_entry(vn, block, vardict, 10 * " ")
                        )
                        s += f"{blockentry}\n"

    return s


def md_replace(s):
    # replace specific latex commands
    re_dict = {
        re.compile("\\\\cite{(.*?)\\}"): ("\\cite{{{}}}", None),
        re.compile("\\\\citep{(.*?)\\}"): ("\\citep{{{}}}", None),
        re.compile("\\\\texttt{(.*?)\\}"): ("\\texttt{{{}}}", "`{}`"),
        re.compile("\\$(.*?)\\$"): ("${}$", "{}"),
        re.compile("\\^{(.*?)\\}"): ("^{{{}}}", "<sup>{}</sup>"),
        re.compile("\\^(.*?)\\ "): ("^{:.1}", "<sup>{:.1}</sup>"),
        re.compile("\\``(.*?)\\''"): ("``{}''", '"{}"'),
        re.compile("\\`(.*?)\\'"): ("`{}'", '"{}"'),
    }
    for key, (in_fmt, out_fmt) in re_dict.items():
        for v in key.findall(s):
            src = in_fmt.format(v)
            if out_fmt is None:
                dst = ""
            else:
                dst = out_fmt.format(v)
            s = s.replace(src, dst)

    # replace individual characters
    replace_dict = {
        "\mf": "MODFLOW 6",
        "~": " ",
        "@": "",
        "$": "",
        "\_": "_",
        "&": "|",
        "\le": "&#8804;",
        "\ge": "&#8805;",
        "\\times": "x",
        "\\": "",
    }
    for key, value in replace_dict.items():
        s = s.replace(key, value)

    # delete lines
    delete_tuple = (r"{tabular}", r"hline")
    for value in delete_tuple:
        if value in s:
            s = ""

    return s


def get_examples(component):
    files = [
        filename
        for filename in sorted(os.listdir(EXAMPLES_DIR_PATH))
        if component.lower() in filename.lower() and "-obs" not in filename.lower()
    ]
    s = ""
    for idx, filename in enumerate(files):
        if idx == 0:
            s += "#### Example Input File\n"
        if len(files) > 1:
            s += f"Example {idx + 1}\n\n"
        fpth = os.path.join(EXAMPLES_DIR_PATH, filename)
        with open(fpth, "r") as f:
            lines = f.readlines()
        s += "```\n"
        for line in lines:
            line = line.rstrip()
            s += f"    {line}\n"
        s += "```\n\n"
    return s


def get_obs_examples(component):
    files = [
        filename
        for filename in sorted(os.listdir(EXAMPLES_DIR_PATH))
        if component.lower() in filename.lower() and "-obs" in filename.lower()
    ]
    s = ""
    for idx, filename in enumerate(files):
        s += "#### Example Observation Input File\n"
        if len(files) > 1:
            s += f"Example {idx + 1}\n\n"
        fpth = os.path.join(EXAMPLES_DIR_PATH, filename)
        with open(fpth, "r") as f:
            lines = f.readlines()
        s += "```\n"
        for line in lines:
            line = line.rstrip()
            s += f"    {line}\n"
        s += "```\n\n"
    return s


def get_obs_table(component):
    files = [
        filename
        for filename in sorted(os.listdir(COMMON_DIR_PATH))
        if component.lower() in filename.lower()
        and filename.lower().endswith("obs.tex")
    ]
    s = ""
    if files:
        s += "#### Available Observation Types\n\n"
        s += "| Stress Package | Observation Type | ID1 | ID2 | Description |\n"
        s += "|----------------|------------------|-----|-----|-------------|\n"
    for filename in files:
        fpth = os.path.join(COMMON_DIR_PATH, filename)
        with open(fpth, "r") as f:
            lines = f.readlines()
        for line in lines:
            line = md_replace(line.rstrip())
            save_line = True
            if len(line) < 1:
                save_line = False
            elif line.strip().startswith("%"):
                save_line = False
            if save_line:
                s += f"| {line.replace('&', '|')} |\n"
    if len(s) > 0:
        s += "\n\n"
    return s


def write_md_header(f):
    s = "# MODFLOW 6 INPUT VARIABLES\n\n"
    f.write(s)
    s = "| component | package | block | variable name | type | description |\n"
    f.write(s)
    s = "| :---: | :---: | :---: | :---: | :---: | --- |\n"
    f.write(s)


def write_md(f, vardict, component, package):
    c = component.upper()
    p = package.upper()
    for name, b in vardict.keys():
        n = name.upper()
        v = vardict[(name, b)]
        b = v["block"].upper()
        t = v["type"].upper()
        s = ""
        if t.startswith("REC"):
            pass
        else:
            if t.startswith("KEYSTRING"):
                t = "KEYSTRING"
            t = f"{t}"
            if "shape" in v:
                shape = v["shape"].upper()
                t = f"{t} {shape}"
            d = get_description(v["description"])
            s = f"| {c} | {p} | {b} | {n} | {t} | {d} |\n"
        f.write(s)


def write_appendix(blocks):
    with open(Path(TEX_DIR_PATH) / "appendixA.tex", "w") as f:
        f.write("\\small\n\\begin{longtable}{p{1.5cm} p{1.5cm} p{3cm} c}\n")
        f.write(
            "\\caption{List of block names organized by component and input file "
            "type.  OPEN/CLOSE indicates whether or not the block information "
            "can be contained in separate file} \\tabularnewline \n\n"
        )
        f.write("\\hline\n\\hline\n")
        f.write(
            "\\textbf{Component} & \\textbf{FTYPE} & \\textbf{Blockname} & "
            "\\textbf{OPEN/CLOSE} \\\\\n"
        )
        f.write("\\hline\n\\endfirsthead\n\n\n")

        f.write("\captionsetup{textformat=simple}\n")
        f.write(
            "\caption*{\\textbf{Table A--\\arabic{table}.}{\quad}List of block"
            " names organized by component and input file type.  OPEN/CLOSE "
            "indicates whether or not the block information can be contained "
            "in separate file.---Continued} \\tabularnewline\n"
        )

        f.write("\n\\hline\n\\hline\n")
        f.write(
            "\\textbf{Component} & \\textbf{FTYPE} & \\textbf{Blockname} & "
            "\\textbf{OPEN/CLOSE} \\\\\n"
        )
        f.write("\\hline\n\\endhead\n\n\\hline\n\\endfoot\n\n\n")

        lastftype = ""
        for b in blocks:
            l = b.strip().split("-")
            component, ftype, blockname = l
            if lastftype != ftype:
                f.write("\\hline\n")
            oc = "yes"
            if "griddata" in blockname.lower():
                oc = "no"
            if (
                "utl" in component.lower()
                and "tas" in ftype.lower()
                and "time" in blockname.lower()
            ):
                oc = "no"
            f.write(
                f"{component.upper()} & {ftype.upper()} & {blockname.upper()} & {oc} "
                "\\\\ \n"
            )
            lastftype = ftype

        f.write("\n\n\\hline\n\\end{longtable}\n\\label{table:blocks}\n\\normalsize\n")


def get_dfn_files(models):
    def is_sim_dfn(stem):
        return "sim" in stem

    def is_sln_dfn(stem):
        return "sln" in stem

    def is_utl_dfn(stem):
        return "utl" in stem

    def is_model_dfn(stem):
        return any(stem.startswith(m) for m in models)

    def is_exg_dfn(stem):
        exg = stem.rpartition("-")[2]
        left = exg[:3]
        right = exg[-3:]
        return left in models and right in models

    files = list(DFNS_DIR_PATH.glob("*.dfn"))
    files = (
        sorted([f for f in files if is_sim_dfn(f.stem)])
        + sorted([f for f in files if is_exg_dfn(f.stem)])
        + sorted([f for f in files if is_sln_dfn(f.stem)])
        + sorted([f for f in files if is_model_dfn(f.stem)])
        + sorted([f for f in files if is_utl_dfn(f.stem)])
    )
    return files


def write_variables():
    allblocks = []  # cumulative list of all block names

    # write markdown input variables file
    with open(MD_DIR_PATH / "mf6ivar.md", "w") as fmd:
        write_md_header(fmd)

        for fpath in dfns:
            component, package = fpath.stem.split("-")[0:2]
            vardict = parse_mf6var_file(fpath)

            # make list of unique block names
            blocks = []
            for k in vardict:
                v = vardict[k]
                b = v["block"]
                if b not in blocks:
                    blocks.append(b)

            # add a full block name to allblocks
            for block in blocks:
                b = f"{component}-{package}-{block}"
                allblocks.append(b)

            # go through each block and write information
            desc = (
                "% DO NOT MODIFY THIS FILE DIRECTLY.  IT IS CREATED BY mf6ivar.py \n\n"
            )
            for b in blocks:
                blk_var_list = []

                # Write the name of the block to the latex file
                desc += f"\\item \\textbf{'{Block: ' + b.upper() + '}'}\n\n"
                desc += "\\begin{description}\n"
                desc += write_desc(vardict, b, blk_var_list, varexcludeprefix="dev_")
                desc += "\\end{description}\n"

                with open(TEX_DIR_PATH / f"{fpath.stem}-{b}.dat", "w") as f:
                    s = (
                        write_block(vardict, b, blk_var_list, varexcludeprefix="dev_")
                        + "\n"
                    )
                    f.write(s)
                    if verbose:
                        print(s)

            with open(TEX_DIR_PATH / f"{fpath.stem}-desc.tex", "w") as f:
                s = desc + "\n"
                f.write(s)
                if verbose:
                    print(s)

            # write markdown description
            mdname = fpath.stem
            with open(RTD_DOC_DIR_PATH / f"{mdname}.md", "w") as f:
                f.write(f"### {mdname.upper()}\n\n")
                f.write("#### Structure of Blocks\n\n")
                f.write("_FOR EACH SIMULATION_\n\n")
                desc = ""
                for b in blocks:
                    blk_var_list = []

                    # Write the name of the block to the latex file
                    desc += f"##### Block: {b.upper()}\n\n"

                    desc += write_desc_md(
                        vardict, b, blk_var_list, varexcludeprefix="dev_"
                    )

                    if "period" in b.lower():
                        f.write("\n_FOR ANY STRESS PERIOD_\n\n")
                    f.write("```\n")
                    s = (
                        md_replace(
                            write_block(
                                vardict,
                                b,
                                blk_var_list,
                                varexcludeprefix="dev_",
                                indent=4,
                            )
                        )
                        + "\n"
                    )
                    f.write(s)
                    f.write("```\n")
                    if verbose:
                        print(s)

                f.write("\n#### Explanation of Variables\n\n")
                f.write(desc)

                # add examples
                s = get_examples(mdname)
                if len(s) > 0:
                    f.write(s)

                # add observation table
                s = get_obs_table(mdname)
                if len(s) > 0:
                    f.write(s)

                # add observation examples
                s = get_obs_examples(mdname)
                if len(s) > 0:
                    f.write(s)

                # special cases
                if "sln-ims" in mdname:
                    with open("../ims_table.tex", "r") as fims:
                        lines = fims.readlines()
                    s = (
                        "\n\n"
                        "#### IMS variable values for the "
                        "available complexity options\n"
                    )
                    for line in lines:
                        line = md_replace(line.rstrip())
                        save_line = True
                        if len(line) < 1:
                            save_line = False
                        elif line.strip().startswith("%"):
                            save_line = False
                        if save_line:
                            if "Variable" in line:
                                prepend = "\n\n"
                                postpend = (
                                    "|----------------"
                                    + "|----------------"
                                    + "|----------------"
                                    + "|----------------|\n"
                                )
                            else:
                                prepend = ""
                                postpend = ""
                            s += f"{prepend}| {line.replace('&', '|')} |\n{postpend}"

                    if len(s) > 0:
                        s += "\n\n"
                        f.write(s)

            # write markdown
            write_md(fmd, vardict, component, package)

    return allblocks


if __name__ == "__main__":
    # parse arguments
    parser = ArgumentParser(
        prog="Generate MF6 IO documentation files from DFN files",
        formatter_class=RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
            Generate markdown and latex files for MF6 IO documents
            from DFN files. This script reads DFN files and creates
            tables and appendices describing input/output variables
            supported by MODFLOW 6.
            """
        ),
    )
    parser.add_argument(
        "-m",
        "--model",
        required=False,
        action="append",
        help="Filter models to include",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        required=False,
        default=False,
        action="store_true",
        help="Whether to show verbose output",
    )
    args = parser.parse_args()
    models = args.model if args.model else DEFAULT_MODELS
    verbose = args.verbose

    # clean/recreate docdir
    if os.path.isdir(RTD_DOC_DIR_PATH):
        shutil.rmtree(RTD_DOC_DIR_PATH)
    os.makedirs(RTD_DOC_DIR_PATH)

    # filter dfn files corresponding to the selected set of models
    # and write variables and appendix to markdown and latex files
    dfns = get_dfn_files(models)
    blocks = write_variables()
    write_appendix(blocks)
    if verbose:
        for block in blocks:
            print(block)
