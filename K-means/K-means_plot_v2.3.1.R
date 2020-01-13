setwd("E:/IMA/Rscript/K-means/")
kmeans_plot<-function(filename="all.deg.genes.fpkm.xls",groupname="group.lst",mincenter=2,maxcenter=5,
                      sortway=1,plotway="all",writelist=T,cellWidth=8, cellHeight=0.01){

  #filename,groupname,mincenter=3,maxcenter=5,sortway=1,plotway="mean"|"integrated"|"all",writelist=T
  filename="all.deg.genes.fpkm.xls"
  groupname="group.lst"
  mincenter=3
  maxcenter=5
  sortway=1
  plotway="all"
  writelist=T
  cellWidth=0.5
  cellHeight=1.5
  library(pheatmap)
  library(RColorBrewer)
  group<-read.delim(groupname,header=T,check.names=F)
  sample_num<-nrow(group)
  sortname<-paste('^',group[,1],sep="",collapse="|")
  data<-read.delim(filename,header = T,check.names=F)
  index <- duplicated(data[,1])
  removedupdata=data[!index,]
  name=removedupdata[,1]
  data = removedupdata[,-1]
  row.names(data)=name
  data<-data[apply(data,1,sum)>1,]
  data<-data[,grep(sortname,names(data))]
  group$num<-0
  for (i in 1:sample_num){
    group[i,]$num<-length(grep(paste0('^',group[i,1]),names(data)))}
  if(0%in%group[,2]){stop("Unrecognized sample names" )}
  rdata<-data
  data<-t(scale(t(data)))
  centers<-c(mincenter:maxcenter)
  #Kmeans
  n<-1
  while(n<=length(centers)){
    ncenter<-centers[n]
    sizelist<-data.frame(matrix(nrow=ncenter,ncol=20))
    clusterlist<-data.frame(matrix(nrow=nrow(data),ncol=20))
    
    uqcluster<-0
    
    while(uqcluster==0){
      for (i in 1:20){
        km<-kmeans(data,center=ncenter,iter.max=ncenter*10)
        sizelist[i]<-sort(km$size)
        clusterlist[i]<-km$cluster
      }
      
      stat_min<-data.frame(table(t(sizelist[1,])))
      min_freq<-max(stat_min[,2])
      minopt_size<-stat_min[stat_min$Freq==min_freq,]
      
      stat_max<-data.frame(table(t(sizelist[ncenter,])))
      max_freq<-max(stat_max[,2])
      maxopt_size<-stat_max[stat_max$Freq==max_freq,]
      
      if(nrow(minopt_size)==1&nrow(maxopt_size)==1){uqcluster<-1}
    }
    
    candidate_cluster<-which(sizelist[1,]==as.numeric(as.character(minopt_size[,1]))&sizelist[ncenter,]==as.numeric(as.character(maxopt_size[,1])))
    #出顺序图
    
    if(length(candidate_cluster)==0){
      n<-n}else{
        
        final_cluster<-clusterlist[,sample(candidate_cluster,1)]
        data_cluster<-cbind(data,final_cluster)
        data_cluster<-data.frame(data_cluster,check.names=F)
        mean_data_cluster<-data.frame(matrix(nrow=nrow(data_cluster),ncol=sample_num))           #求均值表
        for (j in 1:sample_num){
          if(group[j,2]==1){mean_data_cluster[,j]<-data_cluster[,grep(group[j,1],names(data_cluster))]
          }else{mean_data_cluster[,j]<-rowMeans(data_cluster[,grep(group[j,1],names(data_cluster))])}
        }
        if(max(group[,2]==1)&min(group[,2]==1)){plotway="mean"}
        rownames(mean_data_cluster)<-rownames(data_cluster)
        names(mean_data_cluster)<-group[,1]
        mean_data_cluster<-cbind(mean_data_cluster,final_cluster)
        
        #均值填�?
        command<-paste("mean_all_cluster<-","list(",paste(1:ncenter,collapse=","),")",sep="")
        eval(parse(text=command))
        for (k in 1:ncenter){mean_all_cluster[[k]]<-mean_data_cluster[mean_data_cluster$final_cluster==k,]}          #填充all_cluster
        weightlist<-data.frame(matrix(nrow=ncenter,ncol=sample_num))
        names(weightlist)<-group[,1]
        d<-c(sample_num:1,2:sample_num)   #divisor_list
        
        for (a in 1:ncenter){
          for (b in 1:sample_num){
            weightlist[a,b]<-sum(colMeans(mean_all_cluster[[a]][1:sample_num])/d[(sample_num-b+1):(2*sample_num-b)])}   #求权�?
        }		
        sortlist<-data.frame(matrix(nrow=ncenter,ncol=2))
        names(sortlist)<-c("center","sort")
        tmp<-weightlist
        
        if(sortway==1){
          for(s in 1:ncenter){
            sortlist[s,]<-which(weightlist==max(tmp),arr.ind=T)
            tmp_del<-which(tmp==max(tmp),arr.ind=T)
            tmp<-tmp[-tmp_del[1],]
          }}
        
        if(sortway==2){
          for(s in 1:ncenter){
            sortlist[s,]<-which(weightlist==max(tmp),arr.ind=T)
            tmp_del<-which(tmp==max(tmp),arr.ind=T)
            tmp<-tmp[-tmp_del[1],-tmp_del[2]]
          }}	
        sortlist<-sortlist[order(sortlist[,2]),]
        sortlist[,2]<-c(1:ncenter)
        
        ####选择均值或全�?
        road_name<-paste(ncenter,"-centers_k-means",sep="")
        dir.create(road_name)
        
        #均值出�?
        if(plotway=="mean"|plotway=="all"){
          sort_cluster<-mean_data_cluster[1,]
          sort_cluster<-sort_cluster[-1,]
          
          for(f in 1:ncenter){
            sort_cluster<-rbind(sort_cluster,mean_all_cluster[[sortlist[f,1]]])
          }
          #定序
          sort_cluster$sort_cluster<-0
          for(p in 1:ncenter){
            sort_cluster[sort_cluster$final_cluster==sortlist[p,1],]$sort_cluster<-sortlist[p,2]}
          group$num<-1
          plot_cluster<-subset(sort_cluster,select=-c(final_cluster,sort_cluster))
          plot_cluster[plot_cluster>1.5]<-1.5
          plot_cluster[plot_cluster<(-1.5)]<--1.5
          anno_col <- data.frame(GeneClass=factor(rep(group[,1],group[,2])))
          rownames(anno_col) <- colnames(plot_cluster)
          if(sample_num<10){colorlist<-brewer.pal(sample_num,"Set3")}else{colorlist<-rainbow(sample_num)}
          names(colorlist)<-group[,1]
          anno_colors <- list(GeneClass=colorlist)
          gaps_row<-c()
          for(r in 1:ncenter){gaps_row<-c(gaps_row,max(which(sort_cluster$sort_cluster==r)))}
          gaps_col<-c()
          for(c in 1:sample_num){gaps_col<-c(gaps_col,sum(group[1:c,2]))}
          #写图
          #pheatmap(as.matrix(plot_cluster),cluster_row=FALSE,cluster_col=FALSE,fontsize=8,cellwidth=cellWidth*5,cellheight=cellHeight,annotation_col=anno_col,color=colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(100),show_rownames=FALSE,show_colnames=FALSE,annotation_colors=anno_colors,gaps_row=gaps_row,gaps_col=gaps_col,filename=paste(road_name,"/",road_name,"_mean.pdf",sep=""))
          pheatmap(as.matrix(plot_cluster),cluster_row=FALSE,cluster_col=FALSE,fontsize=8,cellwidth=cellWidth*5,cellheight=cellHeight,annotation_col=anno_col,color=colorRampPalette(c("blue","white","red"))(100),show_rownames=FALSE,show_colnames=FALSE,annotation_colors=anno_colors,gaps_row=gaps_row,gaps_col=gaps_col,filename=paste(road_name,"/",road_name,"_mean.pdf",sep=""))

          #pheatmap(as.matrix(plot_cluster),cluster_row=FALSE,cluster_col=FALSE,fontsize=8,width=5.5,height=8,annotation_col=anno_col,color=colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(100),show_rownames=FALSE,show_colnames=FALSE,annotation_colors=anno_colors,gaps_row=gaps_row,gaps_col=gaps_col,filename=paste(road_name,"_mean.pdf",sep=""))
          file.copy(paste(road_name,"/",road_name,"_mean.pdf",sep=""),paste(road_name,"_mean.pdf",sep=""))
        }
        #全值出�?
        if(plotway=="integrated"|plotway=="all"){
          
          sort_cluster<-mean_data_cluster[1,]
          sort_cluster<-sort_cluster[-1,]
          for(f in 1:ncenter){
            sort_cluster<-rbind(sort_cluster,mean_all_cluster[[sortlist[f,1]]])
          }
          #定序
          sort_cluster$sort_cluster<-0
          for(p in 1:ncenter){
            sort_cluster[sort_cluster$final_cluster==sortlist[p,1],]$sort_cluster<-sortlist[p,2]}
          for (i in 1:sample_num){group[i,]$num<-length(grep(paste0('^',group[i,1]),names(rdata)))}
          
          sort_cluster$gene<-rownames(sort_cluster)
          sort_cluster<-sort_cluster[,grep("final_cluster|sort_cluster|gene",names(sort_cluster))]
          data_gene<-data.frame(data,check.names = F)
          data_gene$gene<-rownames(data_gene)
          sort_cluster<-merge(data_gene,sort_cluster,by="gene")
          rownames(sort_cluster)<-sort_cluster$gene
          sort_cluster<-sort_cluster[order(sort_cluster$sort_cluster),]
          sort_cluster<-sort_cluster[,-1]
          plot_cluster<-subset(sort_cluster,select=-c(final_cluster,sort_cluster))
          plot_cluster[plot_cluster>1.5]<-1.5
          plot_cluster[plot_cluster<(-1.5)]<--1.5
          anno_col <- data.frame(GeneClass=factor(rep(group[,1],group[,2])))
          rownames(anno_col)<-colnames(plot_cluster)
          if(sample_num<10){colorlist<-brewer.pal(sample_num,"Set3")}else{colorlist<-rainbow(sample_num)}
          names(colorlist)<-group[,1]
          anno_colors <- list(GeneClass=colorlist)
          gaps_row<-c()
          for(r in 1:ncenter){gaps_row<-c(gaps_row,max(which(sort_cluster$sort_cluster==r)))}
          gaps_col<-c()
          for(c in 1:sample_num){gaps_col<-c(gaps_col,sum(group[1:c,2]))}
          #pheatmap(as.matrix(plot_cluster),cluster_row=FALSE,cluster_col=FALSE,annotation_col=anno_col,color=colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(100),show_rownames=FALSE,show_colnames=FALSE,annotation_colors=anno_colors,gaps_row=gaps_row,gaps_col=gaps_col,filename=paste(road_name,"/",road_name,"_all.pdf",sep=""))
          pheatmap(as.matrix(plot_cluster),cluster_row=FALSE,cluster_col=FALSE,fontsize=8,cellwidth=cellWidth*20,cellheight=cellHeight,annotation_col=anno_col,color=colorRampPalette(c("blue","white","red"))(100),show_rownames=FALSE,show_colnames=FALSE,annotation_colors=anno_colors,gaps_row=gaps_row,gaps_col=gaps_col,filename=paste(road_name,"/",road_name,"_all.pdf",sep=""))
          #pheatmap(as.matrix(plot_cluster),cluster_row=FALSE,cluster_col=FALSE,fontsize=8,width=cellWidth,height=cellHeight,annotation_col=anno_col,color=colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(100),show_rownames=FALSE,show_colnames=FALSE,annotation_colors=anno_colors,gaps_row=gaps_row,gaps_col=gaps_col,filename=paste(road_name,"_all.pdf",sep=""))
          file.copy(paste(road_name,"/",road_name,"_all.pdf",sep=""),paste(road_name,"_all.pdf",sep=""))
        }
        #写表
        if(writelist==T){
          write.table(sort_cluster,paste(road_name,"/",road_name,"_plotdata.txt",sep=""),sep="\t",quote=F,row.names=T)
          gene_cluster<-sort_cluster[,c(ncol(sort_cluster)-1,ncol(sort_cluster))]
          gene_cluster$gene<-rownames(gene_cluster)
          gene_cluster<-gene_cluster[,-1]
          rdata$gene<-rownames(rdata)
          data_out<-merge(rdata,gene_cluster,by="gene")
          data_out<-data_out[order(data_out$sort_cluster),]
          write.table(data_out,paste(road_name,"/",road_name,"_all_cluster.txt",sep=""),sep="\t",quote=F,row.names=F)
          for(w in 1:ncenter){write.table(data_out[data_out$sort_cluster==w,],paste(road_name,"/",road_name,"_cluster-",w,".txt",sep=""),sep="\t",quote=F,row.names=F)}
        }
        n<-n+1
      }
  }
  
}
kmeans_plot()

