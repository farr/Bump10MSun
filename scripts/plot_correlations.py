import numpy as np
import corner

m1pct = np.loadtxt('../paper/figures/m1pct_BPLG__including_230529_3.0.txt')
mu = np.loadtxt('../paper/figures/mu_BPLG__including_230529_3.0.txt')
sigma = np.loadtxt('../paper/figures/sigma_BPLG__including_230529_3.0.txt')

data = np.array([m1pct, mu, sigma]).T

fig = corner.corner(data, labels = [r'$m_{1\%}$', r'$m_{peak}$', r'$\sigma_{peak}$'], plot_datapoints = False, smooth = 1.2)#, levels = [1])
fig.savefig('../paper/figures/correlations.pdf')


