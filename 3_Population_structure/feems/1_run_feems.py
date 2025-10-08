# 1. Import - first import base
import numpy as np
# import pkg_resources -> deprecated
from importlib import resources
from sklearn.impute import SimpleImputer
from pandas_plink import read_plink
import statsmodels.api as sm

# viz
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# feems
from feems.utils import prepare_graph_inputs, cov_to_dist
from feems.objective import comp_mats
from feems.viz import draw_FEEMSmix_surface, plot_FEEMSmix_summary
from feems import SpatialGraph, Objective, Viz

# 2. change matplotlib fonts
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.sans-serif"] = "Arial"

# 2. Read data

data_path = "/home/guillem/Documentos/programari/feems/feems/data/calonectris_100res"

# read the genotype data and mean impute missing data
(bim, fam, G) = read_plink("{}/Calonectris_Cha1pop.chr1.downsampled".format(data_path))
imp = SimpleImputer(missing_values=np.nan, strategy="mean")
genotypes = imp.fit_transform((np.array(G)).T)

# 3. setup graph
coord = np.loadtxt("{}/Calonectris_Cha1pop.chr1.coord".format(data_path))  # sample coordinates
outer = np.loadtxt("{}/Calonectris_Cha1pop.chr1.sh100.outer".format(data_path))  # outer coordinates
#grid_path = "{}/grid_100.shp".format(data_path)  # path to discrete global grid
grid_path = "{}/world_7resolution.shp".format(data_path)  # path to discrete global grid
# graph input files
outer, edges, grid, _ = prepare_graph_inputs(coord=coord, 
                                             ggrid=grid_path,
                                             translated=False, 
                                             buffer=0,
                                             outer=outer)

# 5. Prepare the spatial graph

sp_graph = SpatialGraph(genotypes, coord, grid, edges, scale_snps=True)
#projection = ccrs.AzimuthalEquidistant(central_longitude=5)      # Choose are and projection
# Try and change projection cause Azimutha is crashing

projection=ccrs.Robinson(central_longitude=5, globe=None)

# 6. Exploratory plot - you can skip this part

plt.figure(dpi=200, figsize=(8,6))
#projection = ccrs.AzimuthalEquidistant(central_longitude=-100)
ax = plt.axes(projection=projection)

# add coastlines, borders, and land features
ax.add_feature(cfeature.COASTLINE, edgecolor='#636363', linewidth=0.5)
ax.add_feature(cfeature.BORDERS, edgecolor='gray', linewidth=0.3)
ax.add_feature(cfeature.LAND, facecolor='#f7f7f7')

# add grid points
ax.scatter(grid[:, 0], grid[:, 1], s=3, color='grey', alpha=0.7, transform=ccrs.PlateCarree(), label='grid point')

# add outer boundary
ax.plot(outer[:, 0], outer[:, 1], color='black', linewidth=1, transform=ccrs.PlateCarree(), label='outer boundary')

# add edges
for edge in edges:
    i, j = edge - 1
    ax.plot([grid[i, 0], grid[j, 0]], [grid[i, 1], grid[j, 1]], 
            color='lightgray', linewidth=1, alpha=0.6, transform=ccrs.PlateCarree())

# add sample points
ax.scatter(coord[:, 0], coord[:, 1], s=8, color='black', zorder=2,transform=ccrs.PlateCarree(), label='sample points')
plt.legend()
plt.savefig("prova.png", dpi=200)

# 7. Map of the coordinates and grid

fig = plt.figure(dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=projection)  
v = Viz(ax, sp_graph, projection=projection, edge_width=.5, 
        edge_alpha=1, edge_zorder=100, sample_pt_size=10, 
        obs_node_size=7.5, sample_pt_color="black", 
        cbar_font_size=10)
v.draw_map(longlat=True)
v.draw_samples()
v.draw_edges(use_weights=False)
v.draw_obs_nodes(use_ids=False)
plt.savefig("Calonectris_outer.png")

# 8. Fit feems - but first choose lambda and lambd_q (residual penalization) in the cross validation analysis in script 0_...sh

sp_graph.fit(lamb = 0.001, lamb_q = 0.01, optimize_q = 'n-dim')

# 9. Final plot 

fig = plt.figure(dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=projection)  
v = Viz(ax, sp_graph, projection=projection, edge_width=.5, 
        edge_alpha=1, edge_zorder=100, sample_pt_size=20, 
        obs_node_size=7.5, sample_pt_color="black", 
        cbar_font_size=10)
v.draw_map()
v.draw_edges(use_weights=True)
v.draw_obs_nodes(use_ids=False) 
v.draw_edge_colorbar()

plt.savefig("Calonectris.chr1.downsampled.CONTACT.png")
plt.savefig("Calonectris.chr1.downsampled.CONTACT.pdf")

# 10. Visualize model fit

# creating an obj 
obj = Objective(sp_graph); obj.inv(); obj.grad(reg=False)
# computing distances matrice for fit (expected) vs empirical (observed) 
fit_cov, _, emp_cov = comp_mats(obj)
# subsetting matrices to arrays 
fit_dist = cov_to_dist(fit_cov)[np.tril_indices(sp_graph.n_observed_nodes, k=-1)]
emp_dist = cov_to_dist(emp_cov)[np.tril_indices(sp_graph.n_observed_nodes, k=-1)]

# fitting a linear model to the observed distances
X = sm.add_constant(fit_dist)
mod = sm.OLS(emp_dist, X)
res = mod.fit()
muhat, betahat = res.params

plt.figure(dpi=100)
plt.plot(fit_dist, emp_dist, 'o', color='k', alpha=0.8, markersize=4)
plt.axline((0.5,0.5*betahat+muhat), slope=betahat, color='orange', ls='--', lw=3)
plt.text(1, 0.5, "R²={:.3f}".format(res.rsquared), fontsize=15)
plt.xlabel('Fitted distance'); plt.ylabel('Genetic distance')
plt.title(r"$\tt{FEEMS}$ fit with estimated node-specific variances")

plt.savefig('Calonectris.chr1.downsampled.CONTACT.modelfit.png')
plt.savefig('Calonectris.chr1.downsampled.CONTACT.modelfit.pdf')

