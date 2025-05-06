#https://r-graph-gallery.com/215-the-heatmap-function.html

calo <-read.table("Calonectris.Fst.csv")
heatmap(calo, col = heat.colors(25))

calo_no_Leu <- calo[-7,-7]  # Remove the rows and columns corresponding to Leu
heatmap(calo_no_Leu, col = heat.colors(25))

calo_no_Edw <- calo_no_Leu[-6,-6]  # Remove the rows and columns corresponding to Leu
heatmap(calo_no_Edw, col = heat.colors(50))
