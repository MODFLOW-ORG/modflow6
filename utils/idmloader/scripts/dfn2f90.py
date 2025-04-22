import argparse
from os import PathLike
import sys
import textwrap
from pathlib import Path
from pprint import pprint
from typing import Optional

import yaml
from jinja2 import Environment, FileSystemLoader

from modflow_devtools.dfn import Dfn

from filters import Filters

MF6_LENVARNAME = 16
F90_LINELEN = 82
PROJ_ROOT_PATH = Path(__file__).parents[3]
DEFAULT_DFNS_PATH = Path(__file__).parents[1] / "dfns.txt"
DFN_PATH = PROJ_ROOT_PATH / "doc" / "mf6io" / "mf6ivar" / "dfn"
SRC_PATH = PROJ_ROOT_PATH / "src"
IDM_PATH = SRC_PATH / "Idm"

        

    

    def _write_master(self):
        ofspec = SRC_PATH / "Idm" / "selector" / "IdmDfnSelector.f90"
        with open(ofspec, "w") as fh:
            self._write_master_decl(fh)
            self._write_master_defn(fh, defn="param", dtype="param")
            self._write_master_defn(fh, defn="aggregate", dtype="param")
            self._write_master_defn(fh, defn="block", dtype="block")
            self._write_master_multi(fh)
            self._write_master_sub(fh)
            self._write_master_integration(fh)
            self._write_master_component(fh)
            fh.write("end module IdmDfnSelectorModule\n")

    def _write_selectors(self):
        for c in self._d:
            ofspec = SRC_PATH / "Idm" / "selector" / f"Idm{c.title()}DfnSelector.f90"
            with open(ofspec, "w") as fh:
                self._write_selector_decl(fh, component=c, sc_list=self._d[c])
                self._write_selector_helpers(fh)
                self._write_selector_defn(
                    fh, component=c, sc_list=self._d[c], defn="param", dtype="param"
                )
                self._write_selector_defn(
                    fh, component=c, sc_list=self._d[c], defn="aggregate", dtype="param"
                )
                self._write_selector_defn(
                    fh, component=c, sc_list=self._d[c], defn="block", dtype="block"
                )
                self._write_selector_multi(fh, component=c, sc_list=self._d[c])
                self._write_selector_sub(fh, component=c, sc_list=self._d[c])
                self._write_selector_integration(fh, component=c, sc_list=self._d[c])
                fh.write(f"end module Idm{c.title()}DfnSelectorModule\n")

    def _write_selector_decl(self, fh=None, component=None, sc_list=None):
        space = " "
        c = component
        len_c = len(c)

        s = (
            f"! ** Do Not Modify! MODFLOW 6 system generated file. **\n"
            f"module Idm{c.title()}DfnSelectorModule\n\n"
            f"  use ConstantsModule, only: LENVARNAME\n"
            f"  use SimModule, only: store_error\n"
            f"  use InputDefinitionModule, only: InputParamDefinitionType, &\n"
            f"                                   InputBlockDefinitionType\n"
        )

        for sc in sc_list:
            len_sc = len(sc)
            spacer = space * (len_c + len_sc)

            s += f"  use {c.title()}{sc.title()}InputModule\n"

        s += (
            f"\n  implicit none\n"
            f"  private\n"
            f"  public :: {c.lower()}_param_definitions\n"
            f"  public :: {c.lower()}_aggregate_definitions\n"
            f"  public :: {c.lower()}_block_definitions\n"
            f"  public :: {c.lower()}_idm_multi_package\n"
            f"  public :: {c.lower()}_idm_subpackages\n"
            f"  public :: {c.lower()}_idm_integrated\n\n"
        )
        s += "contains\n\n"

        fh.write(s)

    def _write_selector_helpers(self, fh=None):
        s = (
            "  subroutine set_param_pointer(input_dfn, input_dfn_target)\n"
            "    type(InputParamDefinitionType), dimension(:), "
            "pointer :: input_dfn\n"
            "    type(InputParamDefinitionType), dimension(:), "
            "target :: input_dfn_target\n"
            "    input_dfn => input_dfn_target\n"
            "  end subroutine set_param_pointer\n\n"
        )

        s += (
            "  subroutine set_block_pointer(input_dfn, input_dfn_target)\n"
            "    type(InputBlockDefinitionType), dimension(:), "
            "pointer :: input_dfn\n"
            "    type(InputBlockDefinitionType), dimension(:), "
            "target :: input_dfn_target\n"
            "    input_dfn => input_dfn_target\n"
            "  end subroutine set_block_pointer\n\n"
        )

        s += (
            "  subroutine set_subpkg_pointer(subpkg_list, subpkg_list_target)\n"
            "    character(len=16), dimension(:), "
            "pointer :: subpkg_list\n"
            "    character(len=16), dimension(:), "
            "target :: subpkg_list_target\n"
            "    subpkg_list => subpkg_list_target\n"
            "  end subroutine set_subpkg_pointer\n\n"
        )

        fh.write(s)

    def _write_selector_defn(
        self, fh=None, component=None, sc_list=None, defn=None, dtype=None
    ):
        c = component

        s = (
            f"  function {c.lower()}_{defn.lower()}_definitions(subcomponent) "
            f"result(input_definition)\n"
            f"    character(len=*), intent(in) :: subcomponent\n"
            f"    type(Input{dtype.title()}DefinitionType), dimension(:), "
            f"pointer :: input_definition\n"
            f"    nullify (input_definition)\n"
            f"    select case (subcomponent)\n"
        )

        for sc in sc_list:
            s += (
                f"    case ('{sc}')\n"
                f"      call set_{dtype.lower()}_pointer(input_definition, "
                f"{c.lower()}_{sc.lower()}_{defn.lower()}_definitions)\n"
            )

        s += (
            f"    case default\n"
            f"    end select\n"
            f"    return\n"
            f"  end function {c.lower()}_{defn.lower()}_definitions\n\n"
        )

        fh.write(s)

    def _write_selector_multi(self, fh=None, component=None, sc_list=None):
        c = component

        s = (
            f"  function {c.lower()}_idm_multi_package(subcomponent) "
            f"result(multi_package)\n"
            f"    character(len=*), intent(in) :: subcomponent\n"
            f"    logical :: multi_package\n"
            f"    select case (subcomponent)\n"
        )

        for sc in sc_list:
            s += (
                f"    case ('{sc}')\n"
                f"      multi_package = {c.lower()}_{sc.lower()}_"
                f"multi_package\n"
            )

        s += (
            f"    case default\n"
            f"      call store_error('Idm selector subcomponent "
            f"not found; '//&\n"
            f"                       &'component=\"{c.upper()}\"'//&\n"
            f"                       &', subcomponent=\"'//trim(subcomponent)"
            f"//'\".', .true.)\n"
            f"    end select\n"
            f"    return\n"
            f"  end function {c.lower()}_idm_multi_package\n\n"
        )

        fh.write(s)

    def _write_selector_sub(self, fh=None, component=None, sc_list=None):
        c = component

        s = (
            f"  function {c.lower()}_idm_subpackages(subcomponent) "
            f"result(subpackages)\n"
            f"    character(len=*), intent(in) :: subcomponent\n"
            f"    character(len=16), dimension(:), pointer :: subpackages\n"
            f"    select case (subcomponent)\n"
        )

        for sc in sc_list:
            s += (
                f"    case ('{sc}')\n"
                f"      call set_subpkg_pointer(subpackages, "
                f"{c.lower()}_{sc.lower()}_subpackages)\n"
            )

        s += (
            f"    case default\n"
            f"    end select\n"
            f"    return\n"
            f"  end function {c.lower()}_idm_subpackages\n\n"
        )

        fh.write(s)

    def _write_selector_integration(self, fh=None, component=None, sc_list=None):
        c = component

        s = (
            f"  function {c.lower()}_idm_integrated(subcomponent) "
            f"result(integrated)\n"
            f"    character(len=*), intent(in) :: subcomponent\n"
            f"    logical :: integrated\n"
            f"    integrated = .false.\n"
            f"    select case (subcomponent)\n"
        )

        for sc in sc_list:
            s += f"    case ('{sc}')\n"
            s += "      integrated = .true.\n"

        s += (
            f"    case default\n"
            f"    end select\n"
            f"    return\n"
            f"  end function {c.lower()}_idm_integrated\n\n"
        )

        fh.write(s)

    def _write_master_decl(self, fh=None):
        space = " "

        s = (
            "! ** Do Not Modify! MODFLOW 6 system generated file. **\n"
            "module IdmDfnSelectorModule\n\n"
            "  use ConstantsModule, only: LENVARNAME\n"
            "  use SimModule, only: store_error\n"
            "  use InputDefinitionModule, only: InputParamDefinitionType, &\n"
            "                                   InputBlockDefinitionType\n"
        )

        for c in self._d:
            len_c = len(c)
            spacer = space * (len_c)
            s += f"  use Idm{c.title()}DfnSelectorModule\n"

        s += (
            "\n  implicit none\n"
            "  private\n"
            "  public :: param_definitions\n"
            "  public :: aggregate_definitions\n"
            "  public :: block_definitions\n"
            "  public :: idm_multi_package\n"
            "  public :: idm_subpackages\n"
            "  public :: idm_integrated\n"
            "  public :: idm_component\n\n"
            "contains\n\n"
        )

        fh.write(s)

    def _write_master_defn(self, fh=None, defn=None, dtype=None):
        s = (
            f"  function {defn.lower()}_definitions(component, subcomponent) "
            f"result(input_definition)\n"
            f"    character(len=*), intent(in) :: component\n"
            f"    character(len=*), intent(in) :: subcomponent\n"
            f"    type(Input{dtype.title()}DefinitionType), dimension(:), "
            f"pointer :: input_definition\n"
            f"    nullify (input_definition)\n"
            f"    select case (component)\n"
        )

        for c in dfn_d:
            s += (
                f"    case ('{c}')\n"
                f"      input_definition => {c.lower()}_{defn.lower()}_"
                f"definitions(subcomponent)\n"
            )

        s += (
            f"    case default\n"
            f"    end select\n"
            f"    return\n"
            f"  end function {defn.lower()}_definitions\n\n"
        )

        fh.write(s)

    def _write_master_multi(self, fh=None):
        s = (
            "  function idm_multi_package(component, subcomponent) "
            "result(multi_package)\n"
            "    character(len=*), intent(in) :: component\n"
            "    character(len=*), intent(in) :: subcomponent\n"
            "    logical :: multi_package\n"
            "    select case (component)\n"
        )

        for c in dfn_d:
            s += (
                f"    case ('{c}')\n"
                f"      multi_package = {c.lower()}_idm_multi_"
                f"package(subcomponent)\n"
            )

        s += (
            "    case default\n"
            "      call store_error('Idm selector component not found; '//&\n"
            "                       &'component=\"'//trim(component)//&\n"
            "                       &'\", subcomponent=\"'//trim(subcomponent)"
            "//'\".', .true.)\n"
            "    end select\n"
            "    return\n"
            "  end function idm_multi_package\n\n"
        )

        fh.write(s)

    def _write_master_sub(self, fh=None):
        s = (
            "  function idm_subpackages(component, subcomponent) "
            "result(subpackages)\n"
            "    character(len=*), intent(in) :: component\n"
            "    character(len=*), intent(in) :: subcomponent\n"
            "    character(len=16), dimension(:), pointer :: subpackages\n"
            "    select case (component)\n"
        )

        for c in dfn_d:
            s += (
                f"    case ('{c}')\n"
                f"      subpackages => {c.lower()}_idm_"
                f"subpackages(subcomponent)\n"
            )

        s += (
            "    case default\n"
            "      call store_error('Idm selector component not found; '//&\n"
            "                       &'component=\"'//trim(component)//&\n"
            "                       &'\", subcomponent=\"'//trim(subcomponent)"
            "//'\".', .true.)\n"
            "    end select\n"
            "    return\n"
            "  end function idm_subpackages\n\n"
        )

        fh.write(s)

    def _write_master_integration(self, fh=None):
        s = (
            "  function idm_integrated(component, subcomponent) "
            "result(integrated)\n"
            "    character(len=*), intent(in) :: component\n"
            "    character(len=*), intent(in) :: subcomponent\n"
            "    logical :: integrated\n"
            "    integrated = .false.\n"
            "    select case (component)\n"
        )

        for c in dfn_d:
            s += (
                f"    case ('{c}')\n"
                f"      integrated = {c.lower()}_idm_"
                f"integrated(subcomponent)\n"
            )

        s += (
            "    case default\n"
            "    end select\n"
            "    return\n"
            "  end function idm_integrated\n\n"
        )

        fh.write(s)

    def _write_master_component(self, fh=None):
        s = (
            "  function idm_component(component) "
            "result(integrated)\n"
            "    character(len=*), intent(in) :: component\n"
            "    logical :: integrated\n"
            "    integrated = .false.\n"
            "    select case (component)\n"
        )

        for c in dfn_d:
            s += f"    case ('{c}')\n      integrated = .true.\n"

        s += (
            "    case default\n"
            "    end select\n"
            "    return\n"
            "  end function idm_component\n\n"
        )

        fh.write(s)


def _get_template_env():
    template_loader = FileSystemLoader(Path(__file__).parent)
    template_env = Environment(
    loader=template_loader,
    trim_blocks=True,
    lstrip_blocks=True,
    line_statement_prefix="_",
    keep_trailing_newline=True,
    )
    template_env.filters["value"] = Filters.value
    return template_env


def make_targets(dfn, outdir: PathLike, verbose: bool = False):
    # TODO component file
    # TODO selector file
    pass


def make_all(dfndir: PathLike, outdir: PathLike, verbose: bool = False, version: int = 1):
    """Generate Fortran source files from DFN files."""
    # TODO all component and selector files
    # TODO master selector file
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="Convert DFN files to Fortran source files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
            Generate Fortran source code from DFN files. This script
            converts definition (DFN) files to Fortran source files,
            each representing a parameter set for a particular input
            definition. Fortran files generated by this tool provide
            support for simulations, models or packages described by
            the given DFN files. Each DFN file is transformed into a
            corresponding Fortran file with "idm" and the same stem:
            e.g. gwf-ic.dfn becomes gwf-icidm.f90.
            """
        ),
    )
    parser.add_argument(
        "-d",
        "--dfn",
        required=False,
        default=DEFAULT_DFNS_PATH,
        help="Path to a DFN file, a directory containing DFN files, or to a text or YAML file listing DFN files (one per line)",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        required=False,
        default=IDM_PATH,
        help="The directory to write Fortran source files",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        required=False,
        default=False,
        help="Whether to show verbose output",
    )
    args = parser.parse_args()
    dfn = Path(args.dfn)
    outdir = Path(args.outdir) if args.outdir else Path.cwd()
    verbose = args.verbose

    if dfn.suffix.lower() in [".txt"]:
        dfns = open(dfn, "r").readlines()
        dfns = [l.strip() for l in dfns]
        dfns = [l for l in dfns if not l.startswith("#") and l.lower().endswith(".dfn")]
        if dfn == DEFAULT_DFNS_PATH:
            dfns = [DFN_PATH / p for p in dfns]
    elif dfn.suffix.lower() in [".yml", ".yaml"]:
        dfns = yaml.safe_load(open(dfn, "r"))
    elif dfn.suffix.lower() in [".dfn"]:
        dfns = [dfn]

    assert all(p.is_file() for p in dfns), (
        f"DFNs not found: {[p for p in dfns if not p.is_file()]}"
    )

    if verbose:
        print("Converting DFNs:")
        pprint(dfns)

    

    with open(ofspec, "w") as f:

        param_varnames = [
            var.split(
                f"{self.component.lower()}{self.subcomponent.lower()}_"
            )[1]
            for var in self._param_varnames
        ]

        if not len(self._subpackage):
            self._subpackage.append("".ljust(16))

        f.write(
            self._template.render(
                component=self.component,
                subcomponent=self.subcomponent,
                param_varnames=param_varnames,
                aggregate_varnames=self._aggregate_varnames,
                multi_package=self._multi_package,
                subpackage=self._subpackage,
            )
        )

