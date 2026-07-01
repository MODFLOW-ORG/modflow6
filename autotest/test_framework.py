from framework import TestFramework


def test_run_returns_buffer_when_api_func_fails(monkeypatch, function_tmpdir):
    workspace = function_tmpdir / "libmf6"
    workspace.mkdir()
    target = function_tmpdir / "libmf6.dll"
    target.touch()

    test = TestFramework(
        name="libmf6_failure",
        workspace=function_tmpdir,
        targets={"libmf6": target},
        api_func=lambda *_: (_ for _ in ()).throw(RuntimeError("boom")),
        compare=None,
    )

    monkeypatch.setattr("framework.shutil.which", lambda _: str(target))

    success, buff = test._run(workspace=workspace, target=target)

    assert not success
    assert buff == []
