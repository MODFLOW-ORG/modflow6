#!/usr/bin/env python

"""
This script searches through the mf6 source files and looks for mem_allocate
statements.  These statements are used to build tables that contain
information on the mf6 variables that are stored in the memory manager.

"""

import os

from fortran_parser import get_inheritance_dict, source_dir_to_dict

# Set up and check paths
source_dir = "../../../src"
markdown_file = os.path.join("md", "mf6memvar.md")
latex_file = None
assert os.path.isdir(source_dir)
print("Starting...")


def update_current_class(line, current_class):
    if line.lower().startswith("class") and "this" in line.lower():
        if "(" in line:
            idxstart = line.index("(") + 1
            idxend = line.index(")")
            current_class = line[idxstart:idxend]
            current_class = current_class.strip()
    return current_class


def update_current_module(line, current_module):
    if line.lower().startswith("module"):
        ll = line.strip().split()
        current_module = ll[1]
    return current_module


def parse_mem_allocate_line(line):
    arglist = None
    if "(" in line:
        idxstart = line.index("(") + 1
        idxend = line.rindex(")")
        arglist = line[idxstart:idxend].strip().split(",")
    return arglist


def line_list_to_var_list(line_list, fname):
    "convert clean fortran lines into list of memory manager variables"
    current_class = None
    current_module = None
    class_varname_list = []
    memvar_list = []
    for line in line_list:
        current_class = update_current_class(line, current_class)
        current_module = update_current_module(line, current_module)
        if "mem_allocate(" in line:
            arglist = parse_mem_allocate_line(line)

            # mempath
            memory_path = arglist.pop()

            # variable name
            varname = arglist.pop()
            varname = varname.replace('"', "")
            varname = varname.replace("'", "")

            # fortran variable name
            fortran_varname = arglist.pop(0)
            fortran_varname = fortran_varname.lower()
            npercents = fortran_varname.count("%")
            fortran_varname = fortran_varname.replace("this%", "")
            if current_class is not None:
                class_varname = f"{current_class}.{fortran_varname}"
            else:
                class_varname = fortran_varname

            # source name
            source_name = os.path.basename(fname)

            # shape
            shape = ""
            if len(arglist) > 0:
                shape = "(" + ",".join(arglist) + ")"

            # number of dimensions
            dims = len(arglist)

            # check for uniqueness and write to md and tex
            if class_varname not in class_varname_list:
                class_varname_list.append(class_varname)
                l = [
                    source_name,
                    current_module,
                    current_class,
                    fortran_varname,
                    varname,
                    dims,
                ]
                memvar_list.append(l)

    return memvar_list


def write_md(memvar_list, fmd):
    "write markdown table records for list of memory managed variables"
    for l in memvar_list:
        (source_name, current_module, typename, fortran_varname, varname, dims) = l
        write_md_record(fmd, source_name, current_module, typename, varname, dims)
    return


def write_tex(memvar_list, ftex):
    "write latex table records for list of memory managed variables"
    for l in memvar_list:
        (source_name, current_module, typename, fortran_varname, varname, dims) = l
        write_tex_record(ftex, typename, varname, dims)
    return


def write_md_header(f):
    s = "# MODFLOW 6 MEMORY MANAGER VARIABLES\n\n"
    fmd.write(s)
    s = "| source file | module | type.variable name | variable name | dimensions |\n"
    fmd.write(s)
    s = "| :---: | :---: | :---: | :---: | :---: |\n"
    fmd.write(s)
    return


def write_md_record(f, fname, modulename, classname, varname, varshape):
    s = f"| {fname} | {modulename} | {classname} | {varname} | {varshape} |\n"
    f.write(s)
    return


def write_tex_header(f):
    f.write("\\small\n\\begin{longtable}{p{6cm} p{4cm} p{2cm} }\n")
    f.write(
        "\\caption{List of variables stored in memory manager } \\tabularnewline \n\n"
    )
    f.write("\\hline\n\\hline\n")
    f.write("\\textbf{Class.Variable} & \\textbf{Name} & \\textbf{Dimensions} \\\\\n")
    f.write("\\hline\n\\endfirsthead\n\n\n")

    f.write("\captionsetup{textformat=simple}\n")
    f.write(
        "\caption*{\\textbf{Table B--\\arabic{table}.}{\quad}List of"
        " memory manager variables.---Continued} \\tabularnewline\n"
    )

    f.write("\n\\hline\n\\hline\n")
    f.write("\\textbf{Class.Variable} & \\textbf{Name} & \\textbf{Dimensions} \\\\\n")
    f.write("\\hline\n\\endhead\n\n\\hline\n\\endfoot\n\n\n")


def write_tex_footer(f):
    f.write("\n\n\\hline\n\\end{longtable}\n\\label{table:blocks}\n\\normalsize\n")
    f.close()
    return


def write_tex_record(f, classname, varname, dimension):
    if classname is None:
        classname = ""
    else:
        classname = classname.replace("_", "\_")
        classname = classname.replace("%", "-")
    varname = varname.replace("_", "\_")
    s = f"{classname} & {varname} & {dimension} \\\\ \n"
    f.write(s)
    return


# use fortran parser to return a dictionary with the key equal the file
# name and the value equal to a list of complete fortran lines with comments
# removed and continuation lines appended
d = source_dir_to_dict(source_dir)

# setup a markdown file
fmd = open(markdown_file, "w")
write_md_header(fmd)

# setup a latex file
if latex_file is not None:
    ftex = open(latex_file, "w")
    write_tex_header(ftex)

i = 0
for root, dirs, files in os.walk(source_dir):
    for f in files:
        if f.endswith(".f90"):
            fwpath = os.path.join(root, f)
            full_lines = d[f]
            memvar_list = line_list_to_var_list(full_lines, fwpath)
            if len(memvar_list) > 0:
                print(f"{i} -- {f}")
                i += 1
            write_md(memvar_list, fmd)
            if latex_file is not None:
                write_tex(memvar_list, ftex)

if latex_file is not None:
    write_tex_footer(ftex)

ihd = get_inheritance_dict(source_dir)
print("\ninheritance structure")
for child in ihd:
    print(child, "<--", ihd[child])
