library(ggplot2)

data<-read.csv(snakemake@input[[1]],header=TRUE)

data_fisher<-data[data$Algorithm=="Fisher2011",]

p<-ggplot(data_fisher,aes(factor(N),MCC))+geom_boxplot(aes(fill=factor(degree)))
pdf(snakemake@output[[1]])
print(p)
dev.off()
