from pathlib import Path
import itertools

output_dir = Path('../benchmark')
output_dir.mkdir(exist_ok=True)

densities = [1.e0, 1.e1, 1.e2, 1.e4, 1.e5, 1.e6]
color_temperatures = [5.e3, 1.e4, 2.e4, 4.e4]
luminosities = [1., 10., 100.]

# load cloudy input file template as one big string, and fill in the
# values in curly brackest by treating it as a format string
cloudy_template_file = 'comparison/cloudy/thin-shell-mrn-template/hsphere.in'
with open(cloudy_template_file) as f:
    cloudy_template = f.read();

joblist_file = output_dir / 'joblist.txt'
with open(joblist_file, 'w') as jobf:
    for nh, tc, lum in itertools.product(densities, color_temperatures, luminosities):
        single_point_output_dir = output_dir / '{}_{}_{}'.format(nh, tc, lum)
        single_point_output_dir.mkdir(exist_ok=True)
        cloudy_input = cloudy_template.format(nh=nh, tc=tc, lum=lum)
        with open(single_point_output_dir / 'hsphere.in', 'w') as f:
            f.write(cloudy_input)

        jobf.write('cd {} && cloudy.exe -r hsphere\n'.format(single_point_output_dir))
