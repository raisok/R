library(ggplot2)
library(reshape2)

####
#setwd("e:/IMA/Rscript/DEG_stat/")
# for manual set colour with scale_fill_manual add by yueyao at 2019/11/19
col=c("#E41A1C","#377EB8")
####

result <- dir()[grep("differentially_expressed_genes.xls",dir())]
stat <- matrix(0,length(result),3)
colnames(stat) <- c("Group","Up-regulated","Down-regulated")

for(i in 1:length(result)){
data <- read.csv(result[i],sep="\t",header=T,check.names=F)

stat[i,1] <- sub(".differentially_expressed_genes.xls","",result[i])
stat[i,2] <- as.character(summary(data[,ncol(data)-4] > 0))[3]
stat[i,3] <- as.character(summary(data[,ncol(data)-4] > 0))[2]

mydata <- melt(as.data.frame(stat)[i,],id.vars=1)
mydata$variable <- factor(mydata$variable,levels=c("Up-regulated","Down-regulated"))
mydata$value <- as.numeric(mydata$value)
mydata$Group <- gsub("_vs_"," vs\n",mydata$Group)
out <- ggplot(mydata,aes(x=Group,y=value,fill=variable)) +
  geom_bar(stat="identity",position="dodge",width=0.7) +
  xlab("") +ylab("DEG counts") +labs(fill="Trend") +
  #### modify by yueyao at 2019/11/19
  geom_text(aes(label=value),position=position_dodge(0.7),vjust=0.8,size=(6*5/14)) +
  scale_fill_manual(values=col,limits=c("Up-regulated","Down-regulated"))+theme_classic() +
  #scale_fill_discrete(limits=c("Up-regulated","Down-regulated")) +
  theme_classic() +
  theme(axis.text.x=element_text(hjust=0.5,vjust=0.5,size=6,colour = "black"),
        axis.text.y=element_text(size=6,colour = "black"),
        axis.title=element_text(size=6,colour = "black",face="plain"),
        legend.title=element_text(size=6,colour = "black",face="plain"),
        legend.text = element_text(colour = 'black',face = 'plain',size=6)
        )
#### modify by yueyao at 2019/11/19
ggsave(paste(stat[i,1],".DEG_stat.pdf",sep=""),out,width=4,height=3)
ggsave(paste(stat[i,1],".DEG_stat.png",sep=""),out,width=4,height=3)
}

write.table(stat,file="DEG_stat.xls",row.names=F,quote=F,sep="\t")



mydata <- melt(as.data.frame(stat),id.vars=1)
mydata$variable <- factor(mydata$variable,levels=c("Up-regulated","Down-regulated"))
mydata$value <- as.numeric(mydata$value)
mydata$Group <- gsub("_vs_","vs\n",mydata$Group)
out <- ggplot(mydata,aes(x=Group,y=value,fill=variable)) +
  geom_bar(stat="identity",position="dodge",width=0.7) +
  xlab("") +ylab("DEG counts") +
  labs(fill="Trend") +geom_text(aes(label=value),position=position_dodge(0.7),vjust=0.8,size=(6*5/14)) +
  scale_fill_manual(values=col,limits=c("Up-regulated","Down-regulated"))+theme_classic() +
  #scale_colour_discrete(limits=c("Up-regulated","Down-regulated")) +theme_classic() +
  theme(axis.text.x=element_text(hjust=0.5,vjust=0.5,size=6,colour = "black"),
        axis.text.y=element_text(size=6,colour = "black"),
        axis.title=element_text(size=6,colour = "black",face="plain"),
        legend.title=element_text(size=6,colour = "black",face="plain"),
        legend.text = element_text(colour = 'black',face = 'plain',size=6)
        )
width = nrow(stat) + 2
ggsave(paste("All_comparison.DEG_stat.pdf",sep=""),out,width=width,height=3)
ggsave(paste("All_comparison.DEG_stat.png",sep=""),out,width=width,height=3)
