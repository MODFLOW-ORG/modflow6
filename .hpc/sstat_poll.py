import argparse
import pathlib as pl
import time
from subprocess import PIPE, STDOUT, Popen


def _build_command(
    cmd: str,
    jobid: int,
    fmt: str,
    with_header: bool = True,
) -> list:
    cmd_args = [cmd]
    cmd_args.append(f"-j {jobid}")
    cmd_args.append(f"--format={fmt}")
    if not with_header:
        cmd_args.append("--noheader")
    return cmd_args


def _run_command(
    cmd: str,
    jobid: int,
    end_msg: str,
    fmt: str,
    with_header: bool = True,
    silent: bool = False,
) -> list:
    cmd_args = _build_command(cmd, jobid, fmt, with_header=with_header)
    try:
        result = []
        proc = Popen(cmd_args, stdout=PIPE, stderr=STDOUT, cwd=".")

        while True:
            line = proc.stdout.readline().decode("utf-8")
            if line == "" and proc.poll() is not None:
                break
            line = line.rstrip("\r\n")
            if line:
                if end_msg in line or ".ext+" in line:
                    result = None
                    break
                else:
                    if not silent:
                        print(line)
                    if "-----" not in line:
                        result.append(
                            ",".join(line.split())
                            + ","
                            + time.strftime("%Y-%m-%d %H:%M:%S")
                        )
            else:
                break
    except:
        result = None
    return result


if __name__ == "__main__":
    description = (
        "python script for polling a SLURM job "
        + "while it is running on a fixed interval. "
        + "The python uses the SLURM command 'sstat' "
        + "to return information on the job. By default, "
        + "the script returns JobID, AveCPU, AveRSS, "
        + "and MaxRSS but other data can be returned "
        + "by specifying the format argument "
        + "(--format=JobID,AveCPU,AveRSS,MaxRSS,...)."
    )
    parser = argparse.ArgumentParser("sstat_poll", description=description)
    parser.add_argument("jobid", help="SLURM JobID", type=int)
    parser.add_argument(
        "--format",
        help="SLURM sstat format string (default is JobID,AveCPU,AveRSS,MaxRSS)",
        type=str,
        default="JobID,AveCPU,AveRSS,MaxRSS",
        required=False,
    )
    parser.add_argument(
        "--output",
        help="Output file (default is None)",
        type=str,
        required=False,
        default=None,
    )
    parser.add_argument(
        "--prefix",
        help="Output file prefix (default is None)",
        type=str,
        required=False,
        default=None,
    )
    parser.add_argument(
        "--command",
        help="SLURM function (default is sstat)",
        type=str,
        required=False,
        default="sstat",
    )
    parser.add_argument(
        "--interval",
        help="polling interval in sec. (default is 30.0 sec.) ",
        type=float,
        required=False,
        default=30.0,
    )
    slurm_args = parser.parse_args()

    print(f"SLURM command: {slurm_args.command}")
    print(f"JobID: {slurm_args.jobid}")

    if slurm_args.output is None:
        output_path = f"{slurm_args.jobid}.{slurm_args.command}.csv"
        if slurm_args.prefix is not None:
            output_path = f"{slurm_args.prefix}.{output_path}"
        output_path = pl.Path(output_path)
    else:
        output_path = pl.Path(slurm_args.output)
    print(f"output path: {output_path}")

    end_msg = (
        f"{slurm_args.command}: error: no steps "
        + f"running for job {slurm_args.jobid}"
    )

    # test if exe exists
    if (
        _run_command(
            slurm_args.command,
            slurm_args.jobid,
            end_msg,
            slurm_args.format,
            silent=True,
        )
        is None
    ):
        raise ValueError(f"SLURM command '{slurm_args.command}' does not exist")

    end_tag = f"sstat:,error:,no,steps,running,for,job,{slurm_args.jobid}"
    # open file
    with open(output_path, "w") as f:
        with_header = True
        job_complete = False
        while job_complete is False:
            result = _run_command(
                slurm_args.command,
                slurm_args.jobid,
                end_msg,
                slurm_args.format,
                with_header=with_header,
            )
            if result is None:
                job_complete = True
            if not job_complete:
                with_header = False
                for line in result:
                    f.write(f"{line}\n")
                time.sleep(slurm_args.interval)
