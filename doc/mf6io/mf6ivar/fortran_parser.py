import os


def get_next_line(line_list):
    line = ""
    while len(line_list) > 0:
        line = line_list.pop(0)
        line = rf"{line}"
        line = line.strip()
        if "!" in line:
            idx = line.index("!")
            line = line[:idx]
        line = line.strip()
        if len(line) == 0:
            continue
        else:
            break
    return line, line_list


def get_full_lines(fname):
    with open(fname, "r") as f:
        line_list = f.readlines()
    full_lines = []
    while len(line_list) > 0:
        line, line_list = get_next_line(line_list)
        while line.endswith("&"):
            line = line[:-1].strip()
            lnext, line_list = get_next_line(line_list)
            if lnext.startswith("&"):
                lnext = lnext[1:]
            line += lnext
        full_lines.append(line)
    return full_lines


def source_dir_to_dict(source_dir=".", ext=".f90", verbose=True):
    assert os.path.isdir(source_dir)
    d = {}
    for root, dirs, files in os.walk(source_dir):
        for f in files:
            if f.endswith(ext):
                if verbose:
                    print(f"processing {f}")
                fwpath = os.path.join(root, f)
                full_lines = get_full_lines(fwpath)
                d[f] = full_lines
    return d


def parse_type_extends(line):
    istart = line.index("(")
    istop = line.index(")")
    parent = line[istart + 1 : istop]
    istart = line.index("::")
    child = line[istart + 2 :].strip()
    return child, parent


def get_inheritance_dict(source_dir):
    d = source_dir_to_dict(source_dir)
    inheritance_dict = {}
    for f in d:
        full_lines = d[f]
        for line in full_lines:
            if "type, extends" in line:
                child, parent = parse_type_extends(line)
                inheritance_dict[child] = parent
    return inheritance_dict
