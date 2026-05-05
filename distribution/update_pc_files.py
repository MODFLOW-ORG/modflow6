import pathlib as pl
import subprocess

nc_includepath = pl.Path(
    subprocess.check_output(("nf-config", "--includedir")).decode().strip()
)
nc_libpath = nc_includepath.parent / "lib"
print(
    f"netcdf-fortran include path: {nc_includepath}\n"
    + f"netcdf-fortran lib path:     {nc_libpath}"
)

nc_pcpath = nc_libpath / "pkgconfig/netcdf-fortran.pc"
if nc_pcpath.is_file():
    with open(nc_pcpath, "r") as f:
        lines = f.readlines()
    tag = "fmoddir="
    update = False
    for idx, line in enumerate(lines):
        if tag in line:
            test = line.strip()
            if test == tag:
                update = True
                lines[idx] = line.replace(tag, f"{tag}{nc_includepath}")
            break
    if update:
        print(f"Updating {tag} with the include path.")
        with open(nc_pcpath, "w") as f:
            for line in lines:
                f.write(line)
    else:
        print("All files are in good shape.")
else:
    print(f"netcdf-fortran pkgconfig file '{nc_pcpath}' was not found")
