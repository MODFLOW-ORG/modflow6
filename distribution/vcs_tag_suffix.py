"""
If HEAD is an exact tag match, this is an official release, so
return no suffix. Otherwise it's a development build so return
suffix '+shortsha'.
"""

import subprocess
import sys


def get_suffix():
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
    suffix = get_suffix()
    if len(sys.argv) == 3:
        input_path, output_path = sys.argv[1], sys.argv[2]
        with open(input_path) as f:
            content = f.read().replace("@VCS_TAG@", suffix)
        with open(output_path, "w") as f:
            f.write(content)
    else:
        print(suffix, end="")
