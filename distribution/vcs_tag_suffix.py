"""
Computes the @VCS_TAG@ suffix meson substitutes into version.f90 at
build time. This only ever applies to development builds - anything
not built via `update_version.py --releasemode`, which sets
@VCS_TAG@ to "" itself before meson ever runs.

For a development build, the suffix is '+shortsha[.dirty]', or ''
if HEAD is an exact tag match.
"""

import subprocess


def get_suffix():
    # no vcs tag if git isn't available
    try:
        subprocess.run(["git", "status"], capture_output=True, check=True)
    except (subprocess.CalledProcessError, OSError):
        return ""

    try:
        subprocess.check_output(
            ["git", "describe", "--exact-match", "--tags", "HEAD"],
            stderr=subprocess.DEVNULL,
        )
        return ""
    except subprocess.CalledProcessError:
        sha = (
            subprocess.check_output(["git", "rev-parse", "--short", "HEAD"])
            .decode()
            .strip()
        )
        dirty = (
            subprocess.run(
                ["git", "diff", "--quiet", "HEAD"], capture_output=True
            ).returncode
            != 0
        )
        return f"+{sha}.dirty" if dirty else f"+{sha}"


if __name__ == "__main__":
    print(get_suffix(), end="")
