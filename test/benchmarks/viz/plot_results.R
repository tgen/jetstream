library(ggplot2)
library(dplyr)
library(plotly)
library(reshape2)


df <- data.frame()
results_files <- list.files(pattern = 'trial.*txt')

for(file in results_files){
    data <- read.table(file, sep='', col.names = c('Elapsed', 'CPU', 'Real', 'Virtual'))
    data$Tasks <- as.numeric(gsub('trial[0-9]*_([0-9]{10}).txt', '\\1', file))
    data$Trial <- as.numeric(gsub('trial([0-9]*)_[0-9]{10}.txt', '\\1', file))
    df <- rbind(df, data)
}

plot_data <- df[,c('Trial', 'Tasks', 'Elapsed', 'Real')]

p <- ggplot(plot_data, aes(Elapsed, Real, group=factor(Tasks), color=factor(Tasks)))+
    geom_line(position = 'dodge')
p


plot_data <- df[,c('Trial', 'Tasks', 'Elapsed', 'CPU')]

p <- ggplot(plot_data, aes(Elapsed, CPU, group=factor(Tasks), color=factor(Tasks)))+
    geom_line(position = 'dodge')
p


plot_data <- df[df[,'Tasks'] > 1000, c('Trial', 'Tasks', 'Elapsed', 'Real')]
ggplot(plot_data, aes(Elapsed, Real))+
    geom_point(alpha=0.1)+
    geom_smooth(aes(color=factor(Tasks), fill=factor(Tasks)))


res <- df %>%
        group_by(Tasks, Trial) %>%
        summarise(Elapsed=max(Elapsed),Real=max(Real))

ggplot(res, aes(Tasks, Elapsed, color=factor(Trial))) +
    ggtitle('Elapsed Time vs. Number of Tasks', subtitle = 'Workflow Generation Performance') +
    scale_y_continuous(name='Elapsed Time (s)') +
    geom_line()

ggplot(res, aes(Tasks, Real, color=factor(Trial))) +
    ggtitle('Peak Memory Usage vs. Number of Tasks', subtitle = 'Workflow Generation Performance') +
    scale_y_continuous(name='Real Memory (MB)') +
    geom_line()
