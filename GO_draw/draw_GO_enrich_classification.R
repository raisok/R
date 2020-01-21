setwd("E:\\IMA\\R\\GO_draw")

library(ggplot2)

file="NJFvsLBF_up.DEG_GO_enrichment_result.xls"
outfile="NJFvsLBF_up.DEG_GO_enrichment_result.v2.pdf"

data <- read.table(file,sep="\t",header=T)

data = data[data$Over_represented_pValue <= 0.05,]
data  =  data[data$DEG_item >1,]
#head(data)

shorten_names <- function(x, n_word=4, n_char=30){
  shortname = x["Description"]
  if (length(strsplit(shortname, " ")[[1]]) > n_word || (nchar(shortname) > n_char))
  {
    if (nchar(shortname) > n_char) shortname <- substr(shortname, 1, n_char)
    y <- paste(paste(strsplit(shortname, " ")[[1]][1:min(length(strsplit(shortname," ")[[1]]), n_word)],
                     collapse=" "), "...", sep="")
    return(y)
  } 
  else
  {
    return (shortname)
  }
}

data$pretty_varname = apply(data,1,shorten_names)

draw_data = subset(data,select =c("pretty_varname","Term_type","DEG_item"))

colnames(draw_data) = c("GO_term","GO_category","Num_of_symbols_in_list_in_GO")

draw_data = draw_data[order(draw_data[,2]),]

GO_term_order=factor(draw_data$GO_term,levels=draw_data$GO_term)
COLS=c("red","blue","green")

p <- ggplot(data=draw_data, aes(x=GO_term_order,y=Num_of_symbols_in_list_in_GO, fill=GO_category)) + 
  geom_bar(stat="identity", width=0.8,position=position_dodge(1.0),colour='black',size=0.5) + scale_fill_manual(values=COLS,name="GO Terms",labels=c("BP","CC","MF")) +
  ylab("Number of Genes") + xlab("")+
  labs(title="Enriched GO Terms\n(NJFvsLBF_up)",size=6,face="plain",colour="black")+
  theme_bw()+
  theme(axis.text.x = element_text(color="black",angle=70,hjust = 1,size = 6),
        title =element_text(size=6,face="plain",colour="black"),
        axis.title = element_text(size=6,face="plain",colour="black"),
        panel.grid=element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size=0.2,colour = "black"),
        axis.text = element_text(color = "black",size = 6,face="plain"),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(face="plain",size=6,colour = "black"),
        legend.text = element_text(face="plain",size=6,colour = "black"),
        legend.key.size = unit(0.4,"cm")
        )+scale_x_discrete(expand = c(0, 2))+
  scale_y_discrete(expand = c(0, 0),limits = c(0,5,10,15,20,25,30,35)
                     )
  #guides(fill=FALSE)

ggsave(outfile,p,width=10,height=6)
