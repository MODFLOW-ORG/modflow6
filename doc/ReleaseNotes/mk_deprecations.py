# This script converts the markdown deprecations table
# into a latex table for inclusion in the release notes
import argparse
from pathlib import Path
from warnings import warn

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("path")
    args = parser.parse_args()

    header = r"""
    \section{Deprecations}

    Deprecated/removed options in the current version of MODFLOW 6. Deprecated options are not suggested for use and may (but need not) be removed in a future version of MODFLOW 6. Removed options are no longer available in the current version of MODFLOW 6.

    \small
    \begin{longtable}[!htbp]{p{5cm} p{3cm} p{3cm} p{1.5cm}}
    \caption{List of deprecations and removals}
    \label{table:deprecations}
    \tabularnewline

    \hline
    \hline
    \textbf{Model--Package} & \textbf{Option} & \textbf{Deprecated} & \textbf{Removed}   \\
    \hline
    \endfirsthead

    \hline
    \hline
    \textbf{Model--Package} & \textbf{Option} & \textbf{Deprecated} & \textbf{Removed}   \\
    \hline
    \endhead

    """  # noqa

    footer = r"""
    \hline
    \end{longtable}
    \normalsize
    """

    fname = "deprecations"
    fpath = Path(args.path).expanduser().absolute()
    fnametex = Path(f"{fname}.tex").absolute()
    fnametex.unlink(missing_ok=True)

    # if the markdown file exists, convert it to latex
    if fpath.is_file():
        ftex = open(fnametex, "w")
        ftex.write(header)
        skipline = True
        with open(fpath) as fmd:
            for line in fmd:
                if not skipline:
                    ll = line.strip().split("|")
                    ll = ll[1:-1]
                    linetex = "& ".join(ll)
                    linetex = linetex.replace("\\", "/")
                    linetex += "\\\\" + "\n"
                    linetex = linetex.replace("%", "\\%")
                    linetex = linetex.replace("_", "\\_")
                    ftex.write(linetex)
                    ftex.write("\\hline\n")
                if ":-" in line:
                    skipline = False
        ftex.write(footer)
        ftex.close()
        print(f"Created LaTex file {fnametex} from markdown deprecations file {fpath}")
    else:
        warn(f"Deprecations file not found: {fpath}")
