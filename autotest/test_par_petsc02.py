"""
Test for configuring the PETSc solvers through the petscrc file.
Particularly when running on two different solutions (GWT and GWF)
in parallel.
"""

import os

import pytest
from framework import TestFramework

cases = ["par_petsc02"]


def write_petsc_db(exdir):
    petsc_db_file = os.path.join(exdir, ".petscrc")
    with open(petsc_db_file, "w") as petsc_file:
        petsc_file.write("-SLN_1_ksp_type cg\n")
        petsc_file.write("-SLN_1_use_petsc_pc\n")
        petsc_file.write("-SLN_1_sub_pc_type ilu\n")
        petsc_file.write("-SLN_1_sub_pc_factor_levels 0\n")

        petsc_file.write("-SLN_2_ksp_type bcgs\n")
        petsc_file.write("-SLN_2_use_petsc_pc\n")
        petsc_file.write("-SLN_2_sub_pc_type ilu\n")
        petsc_file.write("-SLN_2_sub_pc_factor_levels 2\n")

    return


def build_models(idx, test):
    from test_gwt_adv01_gwtgwt import build_models as build_models_ext

    # write petscrc
    write_petsc_db(test.workspace)

    return build_models_ext(idx, test)


def check_output(idx, test):
    from test_gwt_adv01_gwtgwt import check_output

    # check model results
    check_output(idx, test)


@pytest.mark.parallel
@pytest.mark.developmode
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        compare=None,
        parallel=True,
        ncpus=2,
    )
    test.run()
