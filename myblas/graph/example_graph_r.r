#import library
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

data <- read.csv(args[1])

plot <- ggplot(data,eas(x=N,y=GFlops))

plot_with_points <- plot + geom_point()

# create png file
png(plot_with_points,filename="plot_with_points.png")

# la fermer/ sauvegarder
dev.off()