library(ggplot2)
library(plotly)
library(reshape2)

df <- data.frame()

for(file in list.files(pattern = 'results.*txt')){
    data <- read.table(file, sep='', col.names = c('Elapsed', 'CPU', 'Real', 'Virtual'))
    data$Tasks <- as.numeric(gsub('results_([0-9]{10}).txt', '\\1', file))
    df <- rbind(df, data)
}

plot_data <- df[,c('Tasks', 'Elapsed', 'Real')]

ggplot(plot_data, aes(Elapsed, Real, color=factor(Tasks)))+
    geom_line()

