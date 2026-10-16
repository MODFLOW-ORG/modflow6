import pytest
from compare import Comparison
from framework import TestFramework

"""
Test for parallel MODFLOW running on two cpus.
It contains two coupled models with

1d:  (nlay,nrow,ncol) = (1,1,5),

constant head boundaries left=1.0, right=10.0.
The result should be a uniform flow field.
"""

cases = ["par_gwf_csv"]
dis_shape = [(1, 1, 5)]

# global convenience...
name_left = "leftmodel"
name_right = "rightmodel"


def update_ims(idx, ims):
    name = cases[idx]
    ims.csv_outer_output_filerecord.set_data(f"{name}.outer.csv")
    ims.csv_inner_output_filerecord.set_data(f"{name}.inner.csv")
    return


def build_models(idx, test):
    from test_par_gwf01 import cases as ex_ext
    from test_par_gwf01 import get_model

    sim = get_model(idx, test.workspace)
    update_ims(idx, sim.get_solution_package(f"{ex_ext[idx]}.ims"))
    return sim, None


def check_output(idx, test):
    from test_par_gwf01 import check_output as check

    check(idx, test)


@pytest.mark.parallel
@pytest.mark.parametrize("idx, name", enumerate(cases))
@pytest.mark.parametrize(
    "compare", [None, pytest.param(Comparison.FP4_ASCII, marks=pytest.mark.fp4)]
)
def test_mf6model(idx, name, compare, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        compare=compare,
        parallel=True,
        ncpus=2,
    )
    test.run()
