# 1. Import
# base
import numpy as np
from importlib import resources
from sklearn.impute import SimpleImputer
from pandas_plink import read_plink
import statsmodels.api as sm

# viz
import matplotlib.pyplot as plt
from matplotlib import gridspec
import cartopy.crs as ccrs

# feems
from feems.utils import prepare_graph_inputs
from feems import SpatialGraph, Viz
from feems.cross_validation import run_cv, run_cv_joint

# change matplotlib fonts
plt.rcParams["font.family"] = "Helvetica"
plt.rcParams["font.sans-serif"] = "Helvetica"

# 2. Read data

data_path = "/home/guillem/Documents/software/feems/feems/data/calonectris"

# read the genotype data and mean impute missing data
(bim, fam, G) = read_plink("{}/Calonectris_Cha1pop.chr1".format(data_path))
imp = SimpleImputer(missing_values=np.nan, strategy="mean")
genotypes = imp.fit_transform((np.array(G)).T)

# 3. setup graph
coord = np.loadtxt("{}/Calonectris_Cha1pop.chr1.coord".format(data_path))  # sample coordinates
outer = np.loadtxt("{}/Calonectris_Cha1pop.chr1.outer".format(data_path))  # outer coordinates
grid_path = "{}/grid_100.shp".format(data_path)  # path to discrete global grid

# graph input files
outer, edges, grid, _ = prepare_graph_inputs(coord=coord, 
                                             ggrid=grid_path,
                                             translated=False, 
                                             buffer=0,
                                             outer=outer)

# 4. construct spatial graph object
sp_graph = SpatialGraph(genotypes, coord, grid, edges, scale_snps=True)

# define grids
# reverse the order of lambdas and alphas for warmstart
lamb_grid = np.geomspace(1e-3, 1e2, 10, endpoint=True)[::-1]
lamb_q_grid = np.geomspace(1e-2, 1e2, 5, endpoint=True)[::-1]

# 5. If you've run the cros_validation_functions.py file, run cross-validation over both smoothing parameters
# ~NEW~ function
cv_err = run_cv_joint(sp_graph, lamb_grid, lamb_q_grid, n_folds=5, factr=1e10)

# average over folds
mean_cv_err = np.mean(cv_err, axis=0)

# lamb & lamb_q CV values
lamb_q_cv = lamb_q_grid[np.where(mean_cv_err == np.min(mean_cv_err))[0][0]]
lamb_cv = lamb_grid[np.where(mean_cv_err == np.min(mean_cv_err))[1][0]]

# plot the figure of CV error
plt.figure(dpi = 200)
plt.gca().set_prop_cycle(color=[plt.get_cmap('Greys_r').resampled(len(lamb_q_grid)+2)(i) for i in range(1,len(lamb_q_grid)+1)])
lineObj = plt.plot(lamb_grid, mean_cv_err.T, '-o', linewidth = 2, alpha = 0.8)
plt.legend(lineObj, lamb_q_grid, title=r'$\lambda_q$'); plt.grid()
plt.xlabel(r'$\lambda$'); plt.semilogx(); plt.ylabel('LOO-CV error')
plt.axvline(lamb_cv, linewidth = 2, color = 'orange')

plt.savefig("prova_7.png")
