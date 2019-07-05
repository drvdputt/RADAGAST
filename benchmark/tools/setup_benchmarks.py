from pathlib import Path
import itertools

c17_dir = Path.home() / 'Software/c17.01/'

gasmodule_git = (Path(__file__).parent / '../..').resolve()

gasmodule_main = (gasmodule_git / '../cmake_release/src/mains/main').resolve()

output_dir = (gasmodule_git / '../benchmark_output').resolve()
output_dir.mkdir(exist_ok=True)

# densities = [1.e0, 1.e1, 1.e2, 1.e3, 1.e4, 1.e5, 1.e6]
densities = [1.e3]
color_temperatures = [1.e4]
luminosities = [10., 100., 1000., 10000., 100000.]

# load cloudy input file template as one big string, and fill in the
# values in curly brackest by treating it as a format string
cloudy_template_file = gasmodule_git / 'benchmark/cloudy/thin-shell-mrn-template/hsphere.in'
with open(cloudy_template_file) as f:
    cloudy_template = f.read();

joblist_file = output_dir / 'joblist.txt'

with open(joblist_file, 'w') as jobf:
    for nh, tc, lum in itertools.product(densities, color_temperatures, luminosities):
        single_point_output_dir = output_dir / '{:.1e}_{:.1e}_{:.1e}'.format(nh, tc, lum)
        single_point_output_dir.mkdir(exist_ok=True)
        cloudy_input = cloudy_template.format(nh=nh, tc=tc, lum=lum)

        with open(single_point_output_dir / 'hsphere.in', 'w') as f:
            f.write(cloudy_input)

        with open(single_point_output_dir / 'parameters.dat', 'w') as f:
            f.writelines((str(d) + '\n' for d in (nh, tc, lum)))

        # job that runs cloudy
        jobf.write('cd {} && {cloudy} -r hsphere\n'.format(single_point_output_dir.resolve(), cloudy=c17_dir / 'source/cloudy.exe'))

        # job that runs gasmodule
        jobf.write('cd {} && mkdir -p MRNDust && {gasmodule} {nh} {tc} {lum}\n'.format(single_point_output_dir.resolve(), gasmodule=gasmodule_main, nh=nh, tc=tc, lum=lum))

run_jobs_command = 'CLOUDY_DATA={}/data parallel < {}'.format(c17_dir, joblist_file)
print('To run the benchmarks, use the following command')
print(run_jobs_command)
