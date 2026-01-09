# Simple script to convert the run time comparison markdown table in the
# distribution folder into a latex table that can be included in the
# release notes
import argparse
from pathlib import Path
from warnings import warn

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("path")
    args = parser.parse_args()

    header = r"""
    \section{Run-Time Comparison}

    Comparison of run times of the current version of MODFLOW 6 to the 
    previous version. The distribution example models are used to compare run 
    times.  Simulations that fail are indicated by '--'. The percent difference, 
    where calculated, is relative to the simulation run time for the previous 
    version. Percent differences for example problems with short run times 
    (less than 30 seconds) may not be significant.

    \small
    \begin{longtable}[!htbp]{p{5cm} p{3cm} p{3cm} p{1.5cm}}
    \caption{List of run time comparsons}
    \label{table:run-time-comparisons}
    \tabularnewline

    \hline
    \hline
    \textbf{Example Problem} & \textbf{Current Version} & \textbf{Previous Version} & \textbf{Percent difference}   \\
    \hline
    \endfirsthead

    \hline
    \hline
    \textbf{Example Problem} & \textbf{Current Version} & \textbf{Previous Version} & \textbf{Percent difference}   \\
    \hline
    \endhead

    """  # noqa

    footer = r"""
    \hline
    \end{longtable}
    \normalsize
    """

    fname = "run-time-comparison"
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
        print(
            f"Created LaTex file {fnametex} "
            f"from markdown benchmark results file {fpath}"
        )
    else:
        warn(f"Benchmark results not found: {fpath}")
