import numpy as np
from pathlib import Path
from matplotlib import pyplot as plt

output_dir = Path('plots')
output_dir.mkdir(exist_ok=True)

# plot grain size distribution
grainabundance = np.loadtxt('grainabundance.g_cm-3.out')
x = range(1, grainabundance.shape[1])
y = grainabundance[0, 1:]
plt.figure()
plt.semilogy(x, y, marker='+')
plt.xlabel('grain pop number')
plt.ylabel('abundance (g_cm-3)')
plt.title('grain abundance')
plt.show()
