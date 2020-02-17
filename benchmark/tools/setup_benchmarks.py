from pathlib import Path
import itertools
import argparse
import glob

gasmodule_git = (Path(__file__).parent / "../..").resolve()

template_dir = gasmodule_git / "benchmark/cloudy"
available_templates = [
    str(Path(s).name).replace("_template", "")
    for s in glob.glob(str(template_dir) + "/*_template")
]

ap = argparse.ArgumentParser()
ap.add_argument(
    "--cloudy-path",
    default=str(Path.home() / "Software/c17.01/"),
    help="default: %(default)s",
)
ap.add_argument(
    "template",
    type=str,
    help="one of the following templates: " + str(available_templates),
)
args = ap.parse_args()

# load cloudy input file template as one big string, and fill in the
# values in curly brackest by treating it as a format string
template_name = args.template
cloudy_template = gasmodule_git / "benchmark/cloudy/{}_template".format(template_name)
if not cloudy_template.is_file():
    raise "{} doesn't exist".format(cloudy_template)

cloudy_input_fname = cloudy_template.name.replace("_template", ".in")

cloudy_dir = Path(args.cloudy_path)
cloudy_exe = cloudy_dir / "source/cloudy.exe"
if not Path(cloudy_exe).is_file():
    raise cloudy_exe + " doesn't exist"

gasmodule_main = (gasmodule_git / "../cmake_release/src/mains/main").resolve()
if not gasmodule_main.is_file():
    raise gasmodule_main + " doesn't exist"

output_dir = (gasmodule_git / "../benchmark_output").resolve()
output_dir.mkdir(exist_ok=True)

densities = [1.0e2, 1.0e3, 1.0e4, 1.0e5]
color_temperatures = [7.5e3, 1.5e4, 3.0e4]
luminosities = [1.0, 10.0, 100.0, 1.0e3, 1.0e4]

with open(cloudy_template) as f:
    cloudy_template_string = f.read()

joblist_file = output_dir / "joblist.txt"
gasmod_jobs_file = output_dir / "only_gasmod.txt"
cloudy_jobs_file = output_dir / "only_cloudy.txt"

with open(joblist_file, "w") as jobf, open(gasmod_jobs_file, "w") as gasmodf, open(
    cloudy_jobs_file, "w"
) as cloudyf:
    for nh, tc, lum in itertools.product(densities, color_temperatures, luminosities):
        output_subdir = output_dir / "{:.1e}_{:.1e}_{:.1e}".format(nh, tc, lum)
        output_subdir.mkdir(exist_ok=True)
        cloudy_input = cloudy_template_string.format(nh=nh, tc=tc, lum=lum)

        with open(output_subdir / cloudy_input_fname, "w") as f:
            f.write(cloudy_input)

        with open(output_subdir / "parameters.dat", "w") as f:
            f.writelines((str(d) + "\n" for d in (nh, tc, lum)))

        # job that runs cloudy
        for f in (jobf, cloudyf):
            f.write(
                "cd {wdir} && {cloudy} -r {name}\n".format(
                    wdir=output_subdir.resolve(), cloudy=cloudy_exe, name=template_name
                )
            )

        # job that runs gasmodule
        for f in (jobf, gasmodf):
            f.write(
                "cd {wdir} && mkdir -p MRNDust && {gasmodule} {nh} {tc} {lum}\n".format(
                    wdir=output_subdir.resolve(),
                    gasmodule=gasmodule_main,
                    nh=nh,
                    tc=tc,
                    lum=lum,
                )
            )

run_jobs_command = "CLOUDY_DATA={}/data parallel < {}".format(cloudy_dir, joblist_file)
print("To run the benchmarks, use the following command")
print(run_jobs_command)
