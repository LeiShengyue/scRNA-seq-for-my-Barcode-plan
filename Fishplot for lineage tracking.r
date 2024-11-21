#load the fishplot library
library(fishplot)

# Define color scheme
mypalette <- c("#D3D3D3", "#b3cde0", "#6497b1", "#005b96", "#011f4b",
               "#8B8378", "#FF00FF", "#8B008B", "#8B5F65", "#BBFFFF", "#00EE76",
               "#FFFACD", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C") # BCD clone colors

# Define timepoints
timepoints = c(1, 2)

# Define parent-child relationships
parents = c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

# Define clone fractions
frac.table = matrix(
    c(
        100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        100, 5.55, 11.11, 5.55, 11.11, 5.55, 5.55, 5.55, 5.55, 5.55, 5.55, 5.55, 5.55, 11.11, 5.55, 5.55),
    ncol = length(timepoints)
)

# Create Fishplot object
fish = createFishObject(frac.table, parents, timepoints = timepoints, fix.missing.clones = TRUE)

# Assign colors to clones
fish2 <- setCol(fish, col = mypalette)

# Layout clones for plotting
fish2 = layoutClones(fish2)

# Generate fishplot (solid background)
pdf("Fishplot.pdf", width = 4, height = 2)
fishPlot(fish2, bg.type = "solid", bg.col = "#FFFFFF", shape = "spline")
dev.off()

# Generate fishplot (gradient background)
pdf("Fishplot2.pdf", width = 4, height = 2)
fishPlot(fish2, bg.type = "solid", bg.col = c("#FFFFFF", "#E0E0E0", "#C0C0C0"), shape = "spline")
dev.off()
