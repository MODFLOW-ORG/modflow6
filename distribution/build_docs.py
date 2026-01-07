import argparse
import os
import platform
import shutil
import sys
import textwrap
from os import PathLike, environ
from pathlib import Path
from pprint import pprint
from tempfile import TemporaryDirectory
from typing import Optional
from urllib.error import HTTPError
from warnings import warn

import pytest
from benchmark import run_benchmarks
from modflow_devtools.build import meson_build
from modflow_devtools.download import (
    download_and_unzip,
    get_release,
)
from modflow_devtools.markers import no_parallel, requires_exe
from modflow_devtools.misc import run_cmd, run_py_script, set_dir

from utils import assert_match, convert_line_endings, get_project_root_path, glob, match

# paths
PROJ_ROOT_PATH = get_project_root_path()
BIN_PATH = PROJ_ROOT_PATH / "bin"
EXAMPLES_REPO_PATH = PROJ_ROOT_PATH.parent / "modflow6-examples"
DISTRIBUTION_PATH = PROJ_ROOT_PATH / "distribution"
DOCS_PATH = PROJ_ROOT_PATH / "doc"
MF6IO_PATH = DOCS_PATH / "mf6io"
MF6IVAR_PATH = MF6IO_PATH / "mf6ivar"
RELEASE_NOTES_PATH = DOCS_PATH / "ReleaseNotes"
TEX_PATHS = {
    "develop": [
        MF6IO_PATH / "mf6io.tex",
        DOCS_PATH / "ReleaseNotes" / "ReleaseNotes.tex",
    ],
    "release": [
        MF6IO_PATH / "mf6io.tex",
        DOCS_PATH / "ReleaseNotes" / "ReleaseNotes.tex",
        DOCS_PATH / "zonebudget" / "zonebudget.tex",
        DOCS_PATH / "ConverterGuide" / "converter_mf5to6.tex",
        DOCS_PATH / "SuppTechInfo" / "mf6suptechinfo.tex",
    ],
}

# models to include in the docs
DEFAULT_MODELS = ["gwf", "gwt", "gwe", "prt"]
DEVELOP_MODELS = ["chf", "olf", "swf"]

# OS-specific extensions
SYSTEM = platform.system()
EXE_EXT = ".exe" if SYSTEM == "Windows" else ""
LIB_EXT = ".dll" if SYSTEM == "Windows" else ".so" if SYSTEM == "Linux" else ".dylib"

# publications
PUB_URLS = [
    "https://pubs.usgs.gov/tm/06/a55/tm6a55.pdf",
    "https://pubs.usgs.gov/tm/06/a56/tm6a56.pdf",
    "https://pubs.usgs.gov/tm/06/a57/tm6a57.pdf",
    "https://pubs.usgs.gov/tm/06/a61/tm6a61.pdf",
    "https://pubs.usgs.gov/tm/06/a62/tm6a62.pdf",
]


@pytest.fixture
def github_user() -> Optional[str]:
    return environ.get("GITHUB_USER", None)


def build_benchmark_tex(
    out_path: PathLike,
    force: bool = False,
    repo_owner: str = "MODFLOW-ORG",
):
    """Build LaTeX files for MF6 performance benchmarks to go into the release notes."""

    # run benchmarks again if no benchmarks found on GitHub or overwrite requested
    benchmarks_path = out_path / "run-time-comparison.md"
    examples_path = EXAMPLES_REPO_PATH / "examples"
    if force or not benchmarks_path.is_file():
        fetch_examples_zip(examples_path, force=force, repo_owner=repo_owner)
        run_benchmarks(
            build_path=PROJ_ROOT_PATH / "builddir",
            current_bin_path=PROJ_ROOT_PATH / "bin",
            previous_bin_path=PROJ_ROOT_PATH / "bin" / "rebuilt",
            examples_path=examples_path,
            out_path=out_path,
        )
    assert benchmarks_path.is_file()

    # convert markdown benchmark results to LaTeX
    with set_dir(RELEASE_NOTES_PATH):
        tex_path = Path("run-time-comparison.tex")
        tex_path.unlink(missing_ok=True)
        out, err, ret = run_cmd(
            sys.executable, "mk_runtimecomp.py", benchmarks_path, verbose=True
        )
        assert not ret, out + err
        assert tex_path.is_file()
    assert (RELEASE_NOTES_PATH / f"{benchmarks_path.stem}.tex").is_file()

    # clean up benchmark results
    benchmarks_path.unlink()


def build_deprecations_tex(force: bool = False):
    """Build LaTeX files for the deprecations table to go into the release notes."""

    # make deprecations markdown table
    (MF6IVAR_PATH / "md").mkdir(exist_ok=True)
    md_path = MF6IVAR_PATH / "md" / "deprecations.md"
    if md_path.is_file() and not force:
        print(f"{md_path} already exists.")
    else:
        md_path.unlink(missing_ok=True)
        with set_dir(MF6IVAR_PATH):
            out, err, ret = run_py_script("deprecations.py", verbose=True)
            assert not ret, out + err

    # convert markdown table to LaTeX
    tex_path = RELEASE_NOTES_PATH / "deprecations.tex"
    if tex_path.is_file() and not force:
        print(f"{tex_path} already exists.")
    else:
        tex_path.unlink(missing_ok=True)
        with set_dir(RELEASE_NOTES_PATH):
            out, err, ret = run_py_script("mk_deprecations.py", md_path, verbose=True)
            assert not ret, out + err

    assert md_path.is_file()
    assert tex_path.is_file()


def build_notes_tex(force: bool = False, patch: bool = False):
    """Build LaTeX files for the release notes."""

    build_deprecations_tex(force=force)

    toml_path = RELEASE_NOTES_PATH / "develop.toml"
    tex_path = RELEASE_NOTES_PATH / "develop.tex"
    if tex_path.is_file() and not force:
        print(f"{tex_path} already exists.")
    else:
        tex_path.unlink(missing_ok=True)
        with set_dir(RELEASE_NOTES_PATH):
            args = [
                "mk_releasenotes.py",
                "--toml",
                toml_path,
                "--tex",
                tex_path,
            ]
            if patch:
                args.append("--patch")
            out, err, ret = run_py_script(*args, verbose=True)
            assert not ret, out + err

    assert tex_path.is_file()


@no_parallel
def test_build_deprecations_tex():
    build_deprecations_tex(force=True)


@no_parallel
def test_build_notes_tex():
    build_notes_tex(force=True)


def build_mf6io_tex(force: bool = False, developmode: bool = True):
    """Build LaTeX files for the MF6IO guide from DFN files."""

    models = DEFAULT_MODELS
    if developmode:
        models.extend(DEVELOP_MODELS)

    included = models + ["sim", "utl", "exg", "sln"]
    excluded = ["appendix", "common"] + list(set(DEFAULT_MODELS) - set(models))

    with set_dir(MF6IVAR_PATH):
        cwd = Path.cwd()

        def _glob(pattern):
            return list(glob(cwd, pattern, included, excluded))

        def _stems(paths):
            return [p.stem.replace("-desc", "") for p in paths]

        tex_files, dfn_files = _glob("*.tex"), _glob("*.dfn")
        tex_stems, dfn_stems = _stems(tex_files), _stems(dfn_files)
        if match(tex_stems, dfn_stems) and not force:
            print("DFN files already exist.")
        else:
            # remove md and tex output dirs
            shutil.rmtree("md", ignore_errors=True)
            shutil.rmtree("tex", ignore_errors=True)

            # run mf6ivar script
            args = [sys.executable, "mf6ivar.py"]
            if not developmode:
                args.append("--releasemode")
            out, err, ret = run_cmd(*args, verbose=True)
            assert not ret, out + err

            # check that a tex file was generated for each dfn
            tex_files, dfn_files = _glob("*.tex"), _glob("*.dfn")
            tex_stems, dfn_stems = _stems(tex_files), _stems(dfn_files)
            assert_match(tex_stems, dfn_stems, "tex", "dfn")


@no_parallel
def test_build_mf6io_tex():
    build_mf6io_tex(force=True)


def build_usage_tex(
    workspace_path: PathLike, bin_path: PathLike, example_model_path: PathLike
):
    """
    Build LaTeX files for the MF6 usage example in the MF6IO guide.
    Runs MF6 to capture the output and insert into the document.
    """

    workspace_path = Path(workspace_path) / "workspace"
    bin_path = Path(bin_path).expanduser().absolute()
    mf6_exe_path = bin_path / f"mf6{EXE_EXT}"
    example_model_path = Path(example_model_path).expanduser().absolute()

    assert mf6_exe_path.is_file(), f"{mf6_exe_path} does not exist"
    assert example_model_path.is_dir(), f"{example_model_path} does not exist"

    tex_path = PROJ_ROOT_PATH / "doc" / "mf6io"
    fname1 = tex_path / "mf6output.tex"
    fname2 = tex_path / "mf6noname.tex"
    fname3 = tex_path / "mf6switches.tex"
    cmd = str(mf6_exe_path)

    if workspace_path.is_dir():
        shutil.rmtree(workspace_path)
    shutil.copytree(example_model_path, workspace_path)

    # run example model
    with set_dir(workspace_path):
        out, err, ret = run_cmd(cmd, verbose=True)
        buff = out + err
        lines = buff.split("\r\n")
        with open(fname1, "w") as f:
            f.write("{\\small\n")
            f.write("\\begin{lstlisting}[style=modeloutput]\n")
            for line in lines:
                f.write(line.rstrip() + "\n")
            f.write("\\end{lstlisting}\n")
            f.write("}\n")

    if workspace_path.is_dir():
        shutil.rmtree(workspace_path)
    os.mkdir(workspace_path)

    # run model without a namefile present
    with set_dir(workspace_path):
        out, err, ret = run_cmd(cmd, verbose=True)
        buff = out + err
        lines = buff.split("\r\n")
        with open(fname2, "w") as f:
            f.write("{\\small\n")
            f.write("\\begin{lstlisting}[style=modeloutput]\n")
            for line in lines:
                f.write(line.rstrip() + "\n")
            f.write("\\end{lstlisting}\n")
            f.write("}\n")

    with set_dir(workspace_path):
        # run mf6 command with -h to show help
        out, err, ret = run_cmd(str(mf6_exe_path), "-h", verbose=True)
        buff = out + err
        lines = buff.split("\r\n")
        with open(fname3, "w") as f:
            f.write("{\\small\n")
            f.write("\\begin{lstlisting}[style=modeloutput]\n")
            for line in lines:
                f.write(line.rstrip() + "\n")
            f.write("\\end{lstlisting}\n")
            f.write("}\n")


def build_pdfs(
    tex_paths: list[PathLike],
    out_path: PathLike,
    passes: int = 3,
    force: bool = False,
):
    """Build PDF documents from LaTeX files."""

    print("Building PDFs from LaTex:")
    pprint(tex_paths)

    out_path = Path(out_path).expanduser().absolute()
    built_paths = set()
    for tex_path in tex_paths:
        tex_path = Path(tex_path).expanduser().absolute()
        pdf_name = tex_path.stem + ".pdf"
        pdf_path = tex_path.parent / pdf_name
        tgt_path = out_path / pdf_name
        if force or not tgt_path.is_file():
            print(f"Converting {tex_path} to PDF")
            with set_dir(tex_path.parent):
                first = True
                for i in range(passes):
                    print(f"Pass {i + 1}/{passes}")
                    out, err, ret = run_cmd(
                        "pdflatex",
                        "-interaction=nonstopmode",
                        "-halt-on-error",
                        tex_path.name,
                    )
                    buff = out + err
                    assert not ret, buff
                    if first:
                        out, err, ret = run_cmd("bibtex", tex_path.stem + ".aux")
                        buff = out + err
                        assert not ret or "I found no" in buff, buff
                        first = False

            if tgt_path.is_file():
                print(f"Clobbering {tgt_path}")
                tgt_path.unlink()

            print(f"Moving {pdf_path} to {tgt_path}")
            pdf_path.rename(tgt_path)
        else:
            print(f"{tgt_path} already exists, nothing to do")

        assert tgt_path.is_file(), f"Failed to build {tgt_path} from {tex_path}"
        assert tgt_path not in built_paths, f"Duplicate target: {tgt_path}"
        built_paths.add(tgt_path)


@no_parallel
@requires_exe("pdflatex")
def test_build_pdfs_from_tex(tmp_path):
    tex_paths = [
        DOCS_PATH / "mf6io" / "mf6io.tex",
        DOCS_PATH / "ReleaseNotes" / "ReleaseNotes.tex",
        DOCS_PATH / "zonebudget" / "zonebudget.tex",
        DOCS_PATH / "ConverterGuide" / "converter_mf5to6.tex",
        DOCS_PATH / "SuppTechInfo" / "mf6suptechinfo.tex",
    ]
    bbl_paths = [
        DOCS_PATH / "ConverterGuide" / "converter_mf5to6.bbl",
    ]

    build_pdfs(tex_paths, tmp_path)

    expected_paths = tex_paths[:-1] + bbl_paths
    assert all(p.is_file() for p in expected_paths)


def fetch_example_docs(
    out_path: PathLike, force: bool = False, repo_owner: str = "MODFLOW-ORG"
):
    pdf_name = "mf6examples.pdf"
    if force or not (out_path / pdf_name).is_file():
        latest = get_release(f"{repo_owner}/modflow6-examples", "latest")
        assets = latest["assets"]
        asset = next(iter([a for a in assets if a["name"] == pdf_name]), None)
        if asset is None:
            raise ValueError(
                f"Release {latest['tag_name']} does not have asset {pdf_name}"
            )
        download_and_unzip(asset["browser_download_url"], out_path, verbose=True)


def fetch_examples_zip(
    out_path: PathLike, force: bool = False, repo_owner: str = "MODFLOW-ORG"
):
    zip_name = "examples.zip"
    out_path = Path(out_path).expanduser().absolute()
    if force or not out_path.is_dir() or not any(os.listdir(out_path)):
        latest = get_release(
            f"{repo_owner}/modflow6-examples", tag="latest", verbose=True
        )
        assets = latest["assets"]
        asset = next(iter([a for a in assets if a["name"].endswith(zip_name)]), None)
        if asset is None:
            raise ValueError(
                f"Release {latest['tag_name']} does not have asset {zip_name}"
            )
        download_and_unzip(asset["browser_download_url"], out_path, verbose=True)


def fetch_usgs_pubs(out_path: PathLike, force: bool = False):
    for url in PUB_URLS:
        print(f"Downloading publication: {url}")
        try:
            download_and_unzip(url, path=out_path, delete_zip=False)
            assert (out_path / url.rpartition("/")[2]).is_file()
        except HTTPError as e:
            if "404" in str(e):
                warn(f"Publication not found: {url}")
            else:
                raise


def build_documentation(
    bin_path: PathLike,
    out_path: PathLike,
    force: bool = False,
    repo_owner: str = "MODFLOW-ORG",
    developmode: bool = True,
    patch: bool = False,
):
    """Build documentation for a MODFLOW 6 distribution."""

    print(f"Building documentation in {'develop' if developmode else 'release'} mode")

    bin_path = Path(bin_path).expanduser().absolute()
    out_path = Path(out_path).expanduser().absolute()
    pdf_path = out_path / "mf6io.pdf"

    if not force and pdf_path.is_file():
        print(f"{pdf_path} already exists, nothing to do")
        return

    out_path.mkdir(parents=True, exist_ok=True)

    with TemporaryDirectory() as temp:
        build_mf6io_tex(force=force, developmode=developmode)
        build_usage_tex(
            bin_path=bin_path,
            workspace_path=Path(temp),
            example_model_path=PROJ_ROOT_PATH / ".mf6minsim",
        )
        build_notes_tex(force=force, patch=patch)

        if developmode:
            tex_paths = TEX_PATHS["develop"]
        else:
            build_benchmark_tex(out_path=out_path, force=force, repo_owner=repo_owner)
            fetch_example_docs(out_path=out_path, force=force, repo_owner=repo_owner)
            fetch_usgs_pubs(out_path=out_path, force=force)
            tex_paths = TEX_PATHS["release"]

        build_pdfs(tex_paths=tex_paths, out_path=out_path, force=force)

    # enforce os line endings on all text files
    windows_line_endings = True
    convert_line_endings(out_path, windows_line_endings)

    # make sure we have expected PDFs
    assert pdf_path.is_file()
    if not developmode:
        assert (out_path / "ReleaseNotes.pdf").is_file()
        assert (out_path / "zonebudget.pdf").is_file()
        assert (out_path / "converter_mf5to6.pdf").is_file()
        assert (out_path / "mf6suptechinfo.pdf").is_file()
        assert (out_path / "mf6examples.pdf").is_file()


@no_parallel
@requires_exe("pdflatex")
def test_build_documentation(tmp_path):
    bin_path = tmp_path / "bin"
    dist_path = tmp_path / "dist"
    meson_build(PROJ_ROOT_PATH, tmp_path / "builddir", bin_path)
    build_documentation(bin_path, dist_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
Create documentation for a distribution. By default, this only includes the mf6io PDF
document. If --releasemode is provided, this includes benchmarks, release notes, the
MODFLOW 6 input/output specification, example model documentation, supplemental info,
documentation for the MODFLOW 5 to 6 converter and Zonebudget 6, and several articles
downloaded from the USGS website too. By default, the script is lazy and will create
only what it can't find. Use the --force (-f) flag to regenerate existing artifacts.
            """
        ),
    )
    parser.add_argument(
        "-b",
        "--bin-path",
        required=False,
        default=str(BIN_PATH),
        help="The path to the directory containing binaries",
    )
    parser.add_argument(
        "-o",
        "--output-path",
        required=False,
        default=os.getcwd(),
        help="The location to create documentation artifacts",
    )
    parser.add_argument(
        "-f",
        "--force",
        required=False,
        default=False,
        action="store_true",
        help="Overwrite existing artifacts. Defaults to false, "
        "so that pre-existing artifacts are used if available.",
    )
    parser.add_argument(
        "--patch",
        default=False,
        action="store_true",
        help="Filter content from release notes for a patch release: "
        "include only items in the 'fixes' section in release notes. "
        "Defaults to false.",
    )
    parser.add_argument(
        "--releasemode",
        required=False,
        default=False,
        action="store_true",
        help="Build documents in release mode for standard releases. "
        "Will omit developmode variables/sections from documentation, "
        "filtering out MF6IO variables marked 'developmode', and also "
        "any LaTeX sections wrapped with '\\ifdevelopmode ... \\fi'. "
        "Defaults false, suitable for preliminary development builds.",
    )
    parser.add_argument(
        "--repo-owner",
        required=False,
        default="MODFLOW-ORG",
        help="Repository owner. Use this option to fetch examples "
        "from a fork of the repository. Defaults to MODFLOW-ORG.",
    )

    args = parser.parse_args()
    bin_path = Path(args.bin_path).expanduser().absolute()
    output_path = Path(args.output_path).expanduser().absolute()
    output_path.mkdir(parents=True, exist_ok=True)
    developmode = not args.releasemode
    repo_owner = args.repo_owner
    force = args.force
    patch = args.patch

    build_documentation(
        bin_path=bin_path,
        out_path=output_path,
        repo_owner=repo_owner,
        developmode=developmode,
        force=force,
        patch=patch,
    )
