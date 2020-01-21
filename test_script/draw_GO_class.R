# 利用下面格式的文件进行作图###############################
# type	terms	number
# CC	plastid	216
# CC	chloroplast	198
# CC	ribosome	151
# CC	extracellular region	83
# CC	thylakoid	60
# CC	photosynthetic membrane	49
# CC	plastid part	98
# CC	thylakoid part	51
# CC	thylakoid membrane	45
# CC	chloroplast part	94
# CC	ribonucleoprotein complex	163
# MF	structural constituent of ribosome	127
# BP	flavonoid metabolic process	13
###########################################################

path = getwd()

setwd("c:/Users/yueyao/Desktop/")
file_in = paste(path,'/go_annotation.txt',sep="")
library(ggplot2)

go <- read.table(file = file_in,sep="\t",header = T)
go$ord=factor(go$terms,levels=go$terms)
mm = 'xxx'
ggplot(data=go,mapping = aes(x=go$ord,y=go$number,fill=go$type))+theme_bw() 
+geom_bar(stat = 'identity',colour='black')+
theme(axis.title = element_text(size=24),axis.title.x =element_text(size=18), 
axis.title.y=element_text(size=18),
axis.text.x = element_text(angle = 30,hjust = 0.98,vjust = 0.98,size=14,))
+labs(x="Number of genes",y="GO terms",title="GO classifcation")+
theme(legend.title = element_blank(),
panel.grid=element_blank(),
legend.position = "top",
legend.text= element_text(size=12),
plot.title = element_text(hjust = 0.5)
)
