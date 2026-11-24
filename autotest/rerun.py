"""
Re-run comparisons: a test simulation is copied, its input re-written, and
the copy run and compared with the original run.

`RERUNS` gives, for each re-run comparison, the steps that re-write and
check a copied simulation workspace. `TestFramework` copies the workspace
(`copy_for_rerun`) and runs MODFLOW 6 between the steps.
"""

import os
import shutil
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from pathlib import Path

import flopy4_input
import mf6_export
import mf6_input
from compare import Comparison, get_mf6_files

# model output, not copied to a re-written simulation's workspace
OUTPUT_EXT = {
    ".lst",
    ".hds",
    ".hed",
    ".bhd",
    ".ahd",
    ".ucn",
    ".cbc",
    ".bud",
    ".grb",
    ".csv",
    ".out",
}


@dataclass(frozen=True)
class Rerun:
    """
    The steps of a re-run comparison, each given a copied simulation
    workspace and, where it takes one, the test's flopy simulation.

    Attributes
    ----------
    prepare : callable
        Re-writes the input, before MODFLOW 6 runs: `prepare(sim_ws, sim)`.
    validate : callable, optional
        If given, MODFLOW 6 first runs in validate mode, then
        `validate(sim_ws, sim)`, before the full run.
    check : callable, optional
        Checks the re-written input or the run's output, after the full
        run: `check(sim_ws)`.
    flopy4 : bool
        Whether the comparison needs flopy4.
    netcdf : bool
        Whether the comparison needs MODFLOW 6 built with NetCDF.
    """

    prepare: Callable
    validate: Callable | None = None
    check: Callable | None = None
    flopy4: bool = False
    netcdf: bool = False


def _flopy4_netcdf(netcdf_format: str) -> Rerun:
    return Rerun(
        prepare=lambda ws, sim: flopy4_input.rewrite_netcdf(ws, netcdf_format),
        check=flopy4_input.check_netcdf_input,
        flopy4=True,
        netcdf=True,
    )


def _mf6_export(netcdf_format: str) -> Rerun:
    return Rerun(
        prepare=lambda ws, sim: mf6_export.enable_export(ws, netcdf_format),
        validate=lambda ws, sim: mf6_export.read_exported(ws, sim, netcdf_format),
        check=mf6_export.check_outputs,
        netcdf=True,
    )


RERUNS = {
    Comparison.FP4_STRUCTURED: _flopy4_netcdf("structured"),
    Comparison.FP4_LAYERED: _flopy4_netcdf("layered"),
    Comparison.FP4_ASCII: Rerun(
        prepare=lambda ws, sim: flopy4_input.rewrite_ascii(ws), flopy4=True
    ),
    Comparison.MF6_EXPORT_STRUCTURED: _mf6_export("structured"),
    Comparison.MF6_EXPORT_LAYERED: _mf6_export("layered"),
}
FLOPY4_COMPARISONS = frozenset(c for c, r in RERUNS.items() if r.flopy4)
NETCDF_COMPARISONS = frozenset(c for c, r in RERUNS.items() if r.netcdf)


def copy_for_rerun(
    src: os.PathLike, dst: os.PathLike, sim_dirs: Iterable[str | os.PathLike]
) -> list[Path]:
    """Copy `src` into `dst` (replaced if it exists) for a re-run, without
    the simulations' output or other comparisons' workspaces. Files a
    re-write doesn't load (e.g. OBS, TS) are kept. Returns `sim_dirs` as
    paths."""
    src, dst = Path(src), Path(dst)
    sim_dir_paths = [Path(d) for d in sim_dirs]
    skip_dirs = {c.value for c in Comparison} - {
        d.parts[0] for d in sim_dir_paths if d.parts
    }

    # output files each simulation names (e.g. a .bin head file)
    named_outputs = {
        d: {Path(f).name for f in get_mf6_files(src / d / "mfsim.nam")[1]}
        for d in sim_dir_paths
    }
    # files an FMI package reads, in any of the simulations (one often
    # reads another's, e.g. "flow/flow.bud")
    fmi_inputs = set().union(
        *(mf6_input.fmi_input_files(src / d) for d in sim_dir_paths)
    )

    def ignore(dirpath, names):
        rel = Path(dirpath).relative_to(src)
        ignored = {n for n in names if rel == Path(".") and n in skip_dirs}
        if rel in sim_dir_paths:
            ignored |= {
                n
                for n in names
                if (Path(n).suffix.lower() in OUTPUT_EXT or n in named_outputs[rel])
                and n not in fmi_inputs
                and (Path(dirpath) / n).is_file()
            }
        return ignored

    if dst.is_dir():
        shutil.rmtree(dst)
    shutil.copytree(src, dst, ignore=ignore)
    return sim_dir_paths
