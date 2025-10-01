# This script doesn't work in its own, it is a continuation of the previous one

# 1. Get the top 10% fraction of outliers (from the distance-genetics matrix?)

outliers_df = sp_graph.extract_outliers(0.01)

# 2. Visualizing the outlier demes on the map
fig = plt.figure(dpi=200)
ax = fig.add_subplot(1, 1, 1, projection=projection)  
v = Viz(ax, sp_graph, projection=projection, edge_width=.5, 
        edge_alpha=1, edge_zorder=100, sample_pt_size=20, 
        obs_node_size=7.5, sample_pt_color="black", 
        cbar_font_size=10)
v.draw_map(); v.draw_edges(use_weights=False)
# ~NEW~ function
v.draw_outliers(outliers_df)
# using deme IDs since all results will be represented with these numbers
v.draw_obs_nodes(use_ids=True)

plt.savefig('prova_5.png')

# 3. Fit long-range edges, from 1 to 10 edges (like in treemix)

seq_results = sp_graph.sequential_fit(
    outliers_df=outliers_df, 
    lamb=2.0, lamb_q=1.0, optimize_q='n-dim', 
    nedges=10, top=5
)      # nedges means k?

# 4. Plot the results.

# How many edges are we plotting here? All?

fig = plt.figure(dpi=250)
ax = fig.add_subplot(1, 1, 1, projection=projection)  
v = Viz(ax, sp_graph, projection=projection, edge_width=.5, 
        edge_alpha=1, edge_zorder=100, sample_pt_size=20, 
        obs_node_size=7.5, sample_pt_color="black", 
        cbar_font_size=10)
v.draw_map(); v.draw_edges(use_weights=True); v.draw_edge_colorbar(); v.draw_obs_nodes()
v.draw_LREs(seq_results); v.draw_c_colorbar()

# 5. Plot and evaluate which is the best k 

plot_FEEMSmix_summary(seq_results, sequential=True)
