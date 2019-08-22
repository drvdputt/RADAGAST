bash -c 'parallel python ../git/benchmark/tools/plot_benchmark_output.py {1}_{2}_1.0e+0? -x lum -prefix n{1}_t{2}_ ::: 1.0e+05 1.0e+04 1.0e+03 1.0e+02 ::: 7.5e+03 1.5e+04 3.0e+04'
