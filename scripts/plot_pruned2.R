# Source: Modified from https://indo-european.eu/human-ancestry/principal-component-analysis-pca-r-eigensoft/
# BEGIN PLOT
# Call to library
#library(rgl)

# select PC components to plot
#xaxe="PC1"
#yaxe="PC2"
args = commandArgs(trailingOnly=TRUE)
xaxe = args[1]
yaxe = args[2]

# Store your xxx.pca.evec file in variable fn
fn <- "out2.evec"

# Read data from fn into data frame evecDat with appropriate column names
evecDat <- read.table(fn, skip = 1, col.names = c("Sample", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "Pop"))

# Read data from study/samples IDs file 
indTable = read.table("st_sample.map", col.names = c("Study","Sample","Symbol","PopColor"))

# Merge both data frames
mergedEvecDat = merge(evecDat, indTable, by = "Sample")

# Now do some corrections : convert to characters colours (if they are text)
mergedEvecDat$PopColor <- as.character(mergedEvecDat$PopColor)

# Read table for studies to use in legend
stTable = read.table("st.dat", col.names = c("Study","Symbol","PopColor"))
stTable$PopColor <- as.character(stTable$PopColor)


# Size of the layout - important for display layout, but also for PDF or images.
layout(matrix(c(1, 2), ncol = 1), heights = c(1.5, 1))
par(mar = c(4, 4, 0, 0))

# open pdf for ploting
pdf(file = paste(xaxe,yaxe,fn,"pdf",sep="."), width = 50,height = 25)
# plot
plot(mergedEvecDat[,eval(xaxe)], mergedEvecDat[,eval(yaxe)],
     main = "PCA using candidate SNPs",
     xlim = c(min(mergedEvecDat[,eval(xaxe)])-0.001, max(mergedEvecDat[,eval(xaxe)])+0.001),
     ylim = c(min(mergedEvecDat[,eval(yaxe)])-0.001, max(mergedEvecDat[,eval(yaxe)])+0.001),
     pch = mergedEvecDat$Symbol,
     col = mergedEvecDat$PopColor,
     cex = 1, cex.axis = 0.9, cex.lab = 0.9,
     xlab = xaxe, ylab = yaxe)

# Write name above samples (uncomment the following)
##text(mergedEvecDat[,eval(xaxe)], mergedEvecDat[,eval(yaxe)], labels = mergedEvecDat$Sample, cex = 0.6, pos = 3)
plot.new()
# Output:
par(mar = rep(0, 4))

# Legend options - this needs correction to show only one name per population, and not the name of all samples:
legend("bottom",
       legend = stTable$Study,
       #other values
       ##legend("center", legend = mergedEvecDat$Study, col = evecDat$Pop, pch = as.integer(evecDat$Pop) %% 24, ncol = 6, cex = 0.6)
       ##legend = levels(indTable$Pop2)
       pch = stTable$Symbol,
       col = stTable$PopColor,
       ncol = 14,
       cex = 0.9)

# Finish plot and print (to PDF, by default)
dev.off()
# END PLOT


