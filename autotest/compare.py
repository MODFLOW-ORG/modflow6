import os
import shutil
from enum import Enum
from pathlib import Path
from typing import Optional, Union
from warnings import warn


class Comparison(Enum):
    """Comparison types for MODFLOW 6 regression testing."""

    AUTO = "auto"
    MF2005 = "mf2005"
    MFNWT = "mfnwt"
    MFUSG = "mfusg"
    MFLGR = "mflgr"
    LIBMF6 = "libmf6"
    MF6 = "mf6"
    MF6_REGRESSION = "mf6_regression"


EXTTEXT = {
    "hds": "head",
    "hed": "head",
    "bhd": "head",
    "ucn": "concentration",
    "cbc": "cell-by-cell",
}
IGNORE_EXTENSIONS = (
    ".hds",
    ".hed",
    ".bud",
    ".cbb",
    ".cbc",
    ".ddn",
    ".ucn",
    ".glo",
    ".lst",
    ".list",
    ".gwv",
    ".mv",
    ".out",
)


def adjust_htol(
    workspace: Union[str, os.PathLike], htol: float = 0.001
) -> Optional[float]:
    """Get outer_dvclose value from MODFLOW 6 ims file"""

    dvclose = get_dvclose(workspace)
    if not dvclose:
        return htol

    # adjust htol if < IMS outer_dvclose
    dvclose *= 5.0
    return dvclose if (htol is None or htol < dvclose) else htol


def get_dvclose(workspace: Union[str, os.PathLike]) -> Optional[float]:
    """Get outer_dvclose value from MODFLOW 6 ims file"""
    dvclose = None
    files = os.listdir(workspace)
    for file_name in files:
        pth = os.path.join(workspace, file_name)
        if os.path.isfile(pth):
            if file_name.lower().endswith(".ims"):
                with open(pth) as f:
                    lines = f.read().splitlines()
                for line in lines:
                    if "outer_dvclose" in line.lower():
                        v = float(line.split()[1])
                        if dvclose is None:
                            dvclose = v
                        else:
                            if v > dvclose:
                                dvclose = v
                        break

    return dvclose


def get_rclose(workspace: Union[str, os.PathLike]) -> Optional[float]:
    """Get inner_rclose value from MODFLOW 6 ims file"""

    rclose = None
    for pth in workspace.glob("*.ims"):
        with open(pth, "r") as f:
            for line in f:
                if "inner_rclose" in line.lower():
                    v = float(line.split()[1])
                    if rclose is None:
                        rclose = v
                    else:
                        if v > rclose:
                            rclose = v
                    break

    if rclose is None:
        return 0.5

    rclose *= 5.0
    return rclose


def get_input_files(namefile):
    """Return a list of all the input files in this model.

    Parameters
    ----------
    namefile : str
        path to a MODFLOW-based model name file

    Returns
    -------
    filelist : list
        list of MODFLOW-based model input files

    """
    srcdir = os.path.dirname(namefile)
    filelist = []
    fname = os.path.join(srcdir, namefile)
    with open(fname, "r") as f:
        lines = f.readlines()

    for line in lines:
        ll = line.strip().split()
        if len(ll) < 2:
            continue
        if line.strip()[0] in ["#", "!"]:
            continue
        ext = os.path.splitext(ll[2])[1]
        if ext.lower() not in IGNORE_EXTENSIONS:
            if len(ll) > 3:
                if "replace" in ll[3].lower():
                    continue
            filelist.append(ll[2])

    # Now go through every file and look for other files to copy,
    # such as 'OPEN/CLOSE'.  If found, then add that file to the
    # list of files to copy.
    otherfiles = []
    for fname in filelist:
        fname = os.path.join(srcdir, fname)
        try:
            f = open(fname, "r")
            for line in f:
                # Skip invalid lines
                ll = line.strip().split()
                if len(ll) < 2:
                    continue
                if line.strip()[0] in ["#", "!"]:
                    continue

                if "OPEN/CLOSE" in line.upper():
                    for i, s in enumerate(ll):
                        if "OPEN/CLOSE" in s.upper():
                            stmp = ll[i + 1]
                            stmp = stmp.replace('"', "")
                            stmp = stmp.replace("'", "")
                            otherfiles.append(stmp)
                            break
        except:
            print(f"{fname} does not exist")

    filelist = filelist + otherfiles

    return filelist


def get_namefiles(pth, exclude=None):
    """Search through a path (pth) for all .nam files.

    Parameters
    ----------
    pth : str
        path to model files
    exclude : str or lst
        File or list of files to exclude from the search (default is None)

    Returns
    -------
    namefiles : lst
        List of namefiles with paths

    """
    namefiles = []
    for root, _, files in os.walk(pth):
        namefiles += [
            os.path.join(root, file) for file in files if file.endswith(".nam")
        ]
    if exclude is not None:
        if isinstance(exclude, str):
            exclude = [exclude]
        exclude = [e.lower() for e in exclude]
        pop_list = []
        for namefile in namefiles:
            for e in exclude:
                if e in namefile.lower():
                    pop_list.append(namefile)
        for e in pop_list:
            namefiles.remove(e)

    return namefiles


def detect_comparison(ws: os.PathLike) -> Comparison | None:
    """
    Determine the comparison type for a MODFLOW 6 regression
    test based on files present in the simulation workspace.
    This function is used to determine which legacy program
    a comparison should be made against, e.g. MODFLOW 2005,
    MODFLOW-NWT, MODFLOW-USG, or MODFLOW-LGR.

    Some files take precedence over others according to the
    order specified in the `Comparison` enum.

    Parameters
    ----------
    src : str
        directory path to search for comparison types

    Returns
    -------
    action : str
        comparison type

    """

    comparisons = [c.value for c in Comparison]
    for _, dirs, _ in os.walk(ws):
        dl = [d.lower() for d in dirs]
        for c in comparisons:
            if any(c in s for s in dl):
                return Comparison(c)
    return None


def get_mf6_files(namefile, verbose=False):
    """Get all MODFLOW 6 input and output files in this simulation.

    Parameters
    ----------
    namefile : pathlike
        path to the MODFLOW 6 simulation name file

    Returns
    -------
    A tuple of lists of paths (input files, output files)
    """

    srcdir = os.path.dirname(namefile)
    mdl_files = []
    pkg_files = []
    out_files = []
    pkg_keys = ["TDIS6", "GWF6", "GWT6", "GWF6-GWF6", "GWF-GWT", "IMS6"]
    model_keys = ["GWF6", "GWT"]

    # find model and simulation-level package input files in simulation namefile
    for line in open(namefile).readlines():
        # Skip over blank and commented lines
        ll = line.strip().split()
        if len(ll) < 2:
            continue
        if line.strip()[0] in ["#", "!"]:
            continue

        for key in pkg_keys:
            if key in ll[0].upper():
                fname = ll[1]
                pkg_files.append(fname)

        for key in model_keys:
            if key in ll[0].upper():
                fname = ll[1]
                mdl_files.append(fname)

    # find model-level package input files in model namefiles
    for namefile in mdl_files:
        fname = os.path.join(srcdir, namefile)
        lines = open(fname, "r").readlines()
        insideblock = False

        for line in lines:
            ll = line.upper().strip().split()
            if len(ll) < 2:
                continue
            if ll[0] in "BEGIN" and ll[1] in "PACKAGES":
                insideblock = True
                continue
            if ll[0] in "END" and ll[1] in "PACKAGES":
                insideblock = False

            if insideblock:
                ll = line.strip().split()
                if len(ll) < 2:
                    continue
                if line.strip()[0] in ["#", "!"]:
                    continue
                pkg_files.append(ll[1])

    # Recurse through package input files and look for input or
    # output file entries, e.g. 'OPEN/CLOSE',  'TIMESERIESFILE'
    # or similar
    flist = pkg_files
    while True:
        olist = []
        flist, olist = get_external_files(srcdir, olist, flist)
        pkg_files += flist
        out_files += olist
        # terminate loop if no additional files
        # if len(flist) < 1 and len(olist) < 1:
        if len(flist) < 1:
            break

    if verbose:
        from pprint import pprint

        print(f"Found input files for {namefile}:")
        pprint(pkg_files)
        print(f"Expecting output files for {namefile}:")
        pprint(out_files)

    return pkg_files, out_files


def get_external_files(srcdir, outplist, files):
    """Get list of external files in a MODFLOW 6 simulation.

    Parameters
    ----------
    srcdir : str
        path to a directory containing a MODFLOW 6 simulation
    outplist : list
        list of output files in a MODFLOW 6 simulation
    files : list
        list of MODFLOW 6 name files

    Returns
    -------

    """
    extfiles = []

    for fname in files:
        fname = os.path.join(srcdir, fname)
        try:
            f = open(fname, "r")
            for line in f:
                # Skip invalid lines
                ll = line.strip().split()
                if len(ll) < 2:
                    continue
                if line.strip()[0] in ["#", "!"]:
                    continue

                if "OPEN/CLOSE" in line.upper():
                    for i, s in enumerate(ll):
                        if s.upper() == "OPEN/CLOSE":
                            stmp = ll[i + 1]
                            stmp = stmp.replace('"', "")
                            stmp = stmp.replace("'", "")
                            extfiles.append(stmp)
                            break

                if "TS6" in line.upper():
                    for i, s in enumerate(ll):
                        if s.upper() == "FILEIN":
                            stmp = ll[i + 1]
                            stmp = stmp.replace('"', "")
                            stmp = stmp.replace("'", "")
                            extfiles.append(stmp)
                            break

                if "TAS6" in line.upper():
                    for i, s in enumerate(ll):
                        if s.upper() == "FILEIN":
                            stmp = ll[i + 1]
                            stmp = stmp.replace('"', "")
                            stmp = stmp.replace("'", "")
                            extfiles.append(stmp)
                            break

                if "OBS6" in line.upper():
                    for i, s in enumerate(ll):
                        if s.upper() == "FILEIN":
                            stmp = ll[i + 1]
                            stmp = stmp.replace('"', "")
                            stmp = stmp.replace("'", "")
                            extfiles.append(stmp)
                            break

                if "EXTERNAL" in line.upper():
                    for i, s in enumerate(ll):
                        if s.upper() == "EXTERNAL":
                            stmp = ll[i + 1]
                            stmp = stmp.replace('"', "")
                            stmp = stmp.replace("'", "")
                            extfiles.append(stmp)
                            break

                if "FILE" in line.upper():
                    for i, s in enumerate(ll):
                        if s.upper() == "FILEIN":
                            stmp = ll[i + 1]
                            stmp = stmp.replace('"', "")
                            stmp = stmp.replace("'", "")
                            extfiles.append(stmp)
                            break

                if "FILE" in line.upper():
                    for i, s in enumerate(ll):
                        if s.upper() == "FILEOUT":
                            stmp = ll[i + 1]
                            stmp = stmp.replace('"', "")
                            stmp = stmp.replace("'", "")
                            outplist.append(stmp)
                            break

        except:
            print("could not get a list of external mf6 files")

    return extfiles, outplist


def get_ftypes(namefile, ftypekeys):
    """Return a list of FTYPES that are in the name file and in ftypekeys.

    Parameters
    ----------
    namefile : str
        path to a MODFLOW 6 name file
    ftypekeys : list
        list of desired FTYPEs

    Returns
    -------
    ftypes : list
        list of FTYPES that match ftypekeys in namefile

    """
    with open(namefile, "r") as f:
        lines = f.readlines()

    ftypes = []
    for line in lines:
        # Skip over blank and commented lines
        ll = line.strip().split()
        if len(ll) < 2:
            continue
        if line.strip()[0] in ["#", "!"]:
            continue

        for key in ftypekeys:
            if ll[0].upper() in key:
                ftypes.append(ll[0])

    return ftypes


def get_comparison_files(
    workspace: os.PathLike, extensions
) -> tuple[list[str], list[str]]:
    if isinstance(extensions, str):
        extensions = [extensions]
    files = os.listdir(workspace)
    files0 = []
    files1 = []
    for file_name in files:
        fpth0 = os.path.join(workspace, file_name)
        if os.path.isfile(fpth0):
            for extension in extensions:
                if file_name.lower().endswith(extension):
                    files0.append(fpth0)
                    fpth1 = os.path.join(workspace, "mf6_regression", file_name)
                    files1.append(fpth1)
                    break
    return files0, files1


def setup_model(namefile, dst, remove_existing=True, extrafiles=None):
    """
    Setup a non-MF6 model test, copying input files to the destination workspace.

    Parameters
    ----------
    namefile : str
        MODFLOW-based model name file.
    dst : str
        destination path for comparison model or file(s)
    remove_existing : bool
        boolean indicating if an existing comparison model or file(s) should
        be replaced (default is True)
    extrafiles : str or list of str
        list of extra files to include in the comparison

    """
    # Construct src pth from namefile or lgr file
    src = os.path.dirname(namefile)

    # Create the destination folder, if required
    create_dir = False
    if os.path.exists(dst):
        if remove_existing:
            print(f"Removing directory '{dst}'")
            shutil.rmtree(dst)
            create_dir = True
    else:
        create_dir = True
    if create_dir:
        os.mkdir(dst)

    # determine if a namefile is a lgr control file - get individual
    # name files out of the lgr control file
    namefiles = [namefile]
    ext = os.path.splitext(namefile)[1]
    if ".lgr" in ext.lower():
        lines = [line.rstrip("\n") for line in open(namefile)]
        for line in lines:
            if len(line) < 1:
                continue
            if line[0] == "#":
                continue
            t = line.split()
            if ".nam" in t[0].lower():
                fpth = os.path.join(src, t[0])
                namefiles.append(fpth)

    # Make list of files to copy
    files2copy = []
    for fpth in namefiles:
        files2copy.append(os.path.basename(fpth))
        ext = os.path.splitext(fpth)[1]
        # copy additional files contained in the name file and
        # associated package files
        if ext.lower() == ".nam":
            fname = os.path.abspath(fpth)
            files2copy = files2copy + get_input_files(fname)

    if extrafiles is not None:
        if isinstance(extrafiles, str):
            extrafiles = [extrafiles]
        for fl in extrafiles:
            files2copy.append(os.path.basename(fl))

    # Copy the files
    for f in files2copy:
        srcf = os.path.join(src, f)
        dstf = os.path.join(dst, f)

        # Check to see if dstf is going into a subfolder, and create that
        # subfolder if it doesn't exist
        sf = os.path.dirname(dstf)
        if not os.path.isdir(sf):
            os.makedirs(sf)

        # Now copy the file
        if os.path.exists(srcf):
            print(f"Copying file '{srcf}' -> '{dstf}'")
            shutil.copy(srcf, dstf)
        else:
            print(f"{srcf} does not exist")


def setup_simulation(
    src, dst, mfnamefile="mfsim.nam", extrafiles=None, remove_existing=True
):
    """
    Setup an MF6 simulation test, copying input files from the source
    to the destination workspace.

    Parameters
    ----------
    src : src
        directory path with original MODFLOW 6 input files
    dst : str
        directory path that original MODFLOW 6 input files will be copied to
    mfnamefile : str
        optional MODFLOW 6 simulation name file (default is mfsim.nam)
    extrafiles : bool
        boolean indicating if extra files should be included (default is None)
    remove_existing : bool
        boolean indicating if existing file in dst should be removed (default
        is True)

    Returns
    -------
    mf6inp : list
        list of MODFLOW 6 input files
    mf6outp : list
        list of MODFLOW 6 output files

    """

    # Create the destination folder
    create_dir = False
    if os.path.exists(dst):
        if remove_existing:
            print(f"Removing {dst}")
            shutil.rmtree(dst)
            create_dir = True
    else:
        create_dir = True
    if create_dir:
        os.makedirs(dst)

    # Make list of files to copy
    fname = os.path.join(src, mfnamefile)
    fname = os.path.abspath(fname)
    mf6inp, mf6outp = get_mf6_files(fname)
    files2copy = [mfnamefile] + mf6inp

    # determine if there are any .ex files
    exinp = []
    for f in mf6outp:
        ext = os.path.splitext(f)[1]
        if ext.lower() == ".hds":
            pth = os.path.join(src, f + ".ex")
            if os.path.isfile(pth):
                exinp.append(f + ".ex")
    if len(exinp) > 0:
        files2copy += exinp
    if extrafiles is not None:
        files2copy += extrafiles

    # Copy the files
    for f in files2copy:
        srcf = os.path.join(src, f)
        dstf = os.path.join(dst, f)

        # Check to see if dstf is going into a subfolder, and create that
        # subfolder if it doesn't exist
        sf = os.path.dirname(dstf)
        if not os.path.isdir(sf):
            try:
                os.mkdir(sf)
            except:
                print(f"Could not create directory '{sf}")

        # Now copy the file
        if os.path.exists(srcf):
            print(f"Copying file '{srcf}' -> '{dstf}'")
            shutil.copy(srcf, dstf)
        else:
            print(f"{srcf} does not exist")

    return mf6inp, mf6outp


def setup_comparison(src, dst, cmp_exe="mf6", overwrite=True, verbose=False):
    """Setup an output comparison for MODFLOW 6 simulation.

    Parameters
    ----------
    src : path-like
        Directory with original MODFLOW 6 input files.
    dst : path-like
        Directory to copy MODFLOW 6 input files to.
    cmp_exe : str or PathLike, optional
        Program to compare with, for supported see the `Comparison` enum.
    overwrite : bool, optional
        Whether to overwrite the destination directory if it exists (default is True).
    verbose : bool, optional
        Whether to show verbose output

    Returns
    -------
    action : str
        comparison type (also the name of the comparison subdirectory in dst)

    """

    if cmp_exe is None:
        warn("No action provided, aborting")
        return

    # create and/or clean dest dir if needed
    dst = Path(dst) / cmp_exe
    dst.mkdir(exist_ok=True)
    dls = list(os.walk(dst))
    if overwrite and any(dls):
        if verbose:
            print(f"Cleaning directory '{dst}'")
        for root, dirs, files in dls:
            for f in files:
                tpth = os.path.join(root, f)
                if verbose:
                    print("Removing file '{tpth}'")
                os.remove(tpth)
            for d in dirs:
                tdir = os.path.join(root, d)
                if verbose:
                    print("Removing directory '{tdir}'")
                shutil.rmtree(tdir)
    else:
        raise ValueError(f"Destination exists but overwrite disabled: {dst}")

    # copy files
    cmppth = os.path.join(src, cmp_exe)
    files = os.listdir(cmppth)
    files2copy = []
    if "mf6" in cmp_exe.lower():
        for file in files:
            if "mfsim.nam" in file.lower():
                srcf = os.path.join(cmppth, os.path.basename(file))
                files2copy.append(srcf)
                srcdir = os.path.join(src, cmp_exe)
                setup_simulation(srcdir, dst, remove_existing=overwrite)
                break
    else:
        for file in files:
            if ".nam" in os.path.splitext(file)[1].lower():
                srcf = os.path.join(cmppth, os.path.basename(file))
                files2copy.append(srcf)
                nf = os.path.join(src, cmp_exe, os.path.basename(file))
                setup_model(nf, dst, remove_existing=overwrite)
                break


def setup_mf5to6(src, dst) -> Path:
    Path(dst).mkdir(exist_ok=True)
    lgrpth = None

    # determine if compare directory exists in directory or if mflgr control
    # file is in directory
    listdir = os.listdir(src)
    for value in listdir:
        fpth = os.path.join(src, value)
        if os.path.isfile(fpth):
            ext = os.path.splitext(fpth)[1]
            if ".lgr" in ext.lower():
                lgrpth = fpth

    print(f"Copying files to target workspace: {dst}")
    # copy lgr files to working directory
    if lgrpth is not None:
        npth = lgrpth
        setup_model(lgrpth, dst)
    # copy MODFLOW-2005, MODFLOW-NWT, or MODFLOW-USG files to working directory
    else:
        npths = get_namefiles(src)
        if len(npths) < 1:
            msg = f"No name files in source workspace: {src}"
            print(msg)
            assert False
        npth = npths[0]
        setup_model(npth, dst)

    return Path(npth)
