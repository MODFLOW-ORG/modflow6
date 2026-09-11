"""
No suffix for release builds. For development builds, suffix '+shortsha[.dirty]'.
"""

import re
import subprocess
import sys


def is_release_build(input_path):
    with open(input_path) as f:
        content = f.read()
    match = re.search(r"IDEVELOPMODE\s*=\s*(\d)", content)
    return match is not None and match.group(1) == "0"


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
    if len(sys.argv) not in (1, 3):
        print(f"usage: {sys.argv[0]} [input output]", file=sys.stderr)
        sys.exit(1)
    if len(sys.argv) == 3:
        input_path, output_path = sys.argv[1], sys.argv[2]
        suffix = "" if is_release_build(input_path) else get_suffix()
        with open(input_path) as f:
            content = f.read().replace("@VCS_TAG@", suffix)
        with open(output_path, "w") as f:
            f.write(content)
    else:
        print(get_suffix(), end="")
