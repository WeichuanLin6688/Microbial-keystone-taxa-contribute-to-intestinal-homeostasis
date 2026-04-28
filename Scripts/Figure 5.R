#fig. 5 ----
#fig. 5B ----
##AB-PAS count-----
df <- read.csv("~/Desktop/Figure 5/CS ABS-PAS count.csv",check.names = F)
df$group <- factor(df$group,levels = c('Young','Old','Old+C. s'))
c <- with(df,kruskal(value,group,p.adj = 'BH'))
c1 <- c$groups
c2 <- c$means
c3 <- merge(c1[,-1,drop=F],c2,by='row.names')
colnames(c3)[1] <- 'group'
c3$group <- factor(c3$group,levels = c('Young','Old','Old+C. s'))
pabpas <- ggplot(df,aes(x=group,y=value,fill=group))+
	geom_point(shape=21,size=2,alpha = 1,stroke=0.3,color='black',
						 position=position_jitterdodge(0.8),show.legend = F)+
	geom_bar(data=c3,stat="identity",width = 0.60,fill='transparent',color='black',lwd=0.2)+
	geom_errorbar(data=c3,aes(x=group,y=value,ymax = value+std,ymin = value-std),lwd=0.2,
								position = position_dodge(0.8),
								width = 0.3,color='black')+
	scale_fill_manual(values = c('transparent','#E24349','#3E7EB2'))+
	labs(y="Cell counts\nper millimeter crypt",x=NULL,title = 'Goblet cells')+
	theme_prism(base_size = 8,	base_family = "sans",
							base_fontface = "plain",base_line_size = 0.2,
							base_rect_size = 0.3,axis_text_angle = 0,border = F )+
	theme(plot.title = element_text(size=9,face ='plain',color='black',margin = margin(b=5)),
				legend.position = 'right',
				axis.text.x = element_text(size=8,angle = 45,hjust = 1,vjust = 1),
				panel.background = element_blank(),plot.background = element_blank())		+
	geom_signif(comparisons = list(c('Young','Old'),c('Old','Old+C. s')),
							map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
							y_position = c(45,50),
							tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
							test = "t.test",test.args = 'two.sided')+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,60),breaks = seq(0,60,20))
pabpas 
#ggsave('~/Desktop/pCS AB-PAS count.pdf',pabpas,width = 1.4,height = 1.80,dpi = 600,units = 'in')

##γH2AX----
df <- subset(read.csv("~/Desktop/Figure 5/H2AX and P16.csv",check.names = F),!name=='P16')
df$group <- factor(df$group,levels = c('Young','Old','Old+C.s'))
c <- with(df,kruskal(mean,group,p.adj = 'BH'))
c1 <- c$groups
c2 <- c$means
c3 <- merge(c1[,-1,drop=F],c2,by='row.names')
colnames(c3)[1] <- 'group'
c3$group <- factor(c3$group,levels = c('Young','Old','Old+C.s'))
df$group <- factor(df$group,levels = c('Young','Old','Old+C.s'))

pH2AX <- ggplot(df,aes(x=group,y=mean,fill=group))+
	geom_point(shape=21,size=2,alpha = 1,stroke=0.3,color='black',
						 position=position_jitterdodge(0.8),show.legend = T)+
	geom_bar(data=c3,stat="identity",width = 0.6,fill='transparent',color='black',lwd=0.2)+
	geom_errorbar(data=c3,aes(x=group,y=mean,ymax = mean+std,ymin = mean-std),lwd=0.2,
								position = position_dodge(0.8),
								width = 0.3,color='black')+
	scale_fill_manual(values = c('transparent','#E24349','#3E7EB2'))+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,2.4),breaks = seq(0,2.4,1.2))+
	labs(y="Area (%)",x=NULL,title = 'γH2AX')+
	theme_prism(base_size = 8,	base_family = "sans",
							base_fontface = "plain",base_line_size = 0.2,
							base_rect_size = 0.3,axis_text_angle = 0,border = F )+
	theme(plot.title = element_text(size=8,face ='plain',color='black',margin = margin(b=5)),
				legend.position = 'none',
				axis.text.x = element_text(size=8,angle = 45,hjust = 1,vjust = 1),
				panel.background = element_blank(),plot.background = element_blank())+
	geom_signif(data=df,aes(x=group,y=mean),
							comparisons = list(c('Young','Old'),c('Old','Old+C.s')),
							map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
							y_position = c(1.9,2.05),
							tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
							test = "t.test",test.args = 'two.sided')
pH2AX
#ggsave('~/Desktop/pH2AX.pdf',pH2AX,width = 1.35,height = 1.8,dpi = 600,units = 'in')

##P16----
df <- subset(read.csv("~/Desktop/Figure 5/H2AX and P16.csv",check.names = F),name=='P16')
df$group <- factor(df$group,levels = c('Young','Old','Old+C.s'))
c <- with(df,kruskal(mean,group,p.adj = 'BH'))
c1 <- c$groups
c2 <- c$means
c3 <- merge(c1[,-1,drop=F],c2,by='row.names')
colnames(c3)[1] <- 'group'
c3$group <- factor(c3$group,levels = c('Young','Old','Old+C.s'))
df$group <- factor(df$group,levels = c('Young','Old','Old+C.s'))

p16 <- ggplot(df,aes(x=group,y=mean,fill=group))+
	geom_point(shape=21,size=2,alpha = 1,stroke=0.3,color='black',
						 position=position_jitterdodge(0.8),show.legend = T)+
	geom_bar(data=c3,stat="identity",width = 0.6,fill='transparent',color='black',lwd=0.2)+
	geom_errorbar(data=c3,aes(x=group,y=mean,ymax = mean+std,ymin = mean-std),lwd=0.2,
								position = position_dodge(0.8),
								width = 0.3,color='black')+
	scale_fill_manual(values = c('transparent','#E24349','#3E7EB2'))+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,2),breaks = seq(0,2,1))+
	labs(y="Area (%)",x=NULL,title = 'P16')+
	theme_prism(base_size = 8,	base_family = "sans",
							base_fontface = "plain",base_line_size = 0.2,
							base_rect_size = 0.3,axis_text_angle = 0,border = F )+
	theme(plot.title = element_text(size=8,face ='plain',color='black',margin = margin(b=5)),
				legend.position = 'none',
				axis.text.x = element_text(size=8,angle = 45,hjust = 1,vjust = 1),
				panel.background = element_blank(),plot.background = element_blank())+
	geom_signif(data=df,aes(x=group,y=mean),
							comparisons = list(c('Young','Old'),c('Old','Old+C.s')),
							map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
							y_position = c(1.6,1.7),
							tip_length = 0,color='black',textsize = 3,vjust = -0.2,size = 0.2,
							test = "t.test",test.args = 'two.sided')
p16
#ggsave('~/Desktop/p16.pdf',p16,width = 1.35,height = 1.8,dpi = 600,units = 'in')


#fig. 5C ----
##aged qPCR------
df1 <- read.csv("~/Desktop/Figure 5/CS aged qPCR.csv",check.names = F)
n <- unique(df1$name);
res <- NULL
for (i in 1:length(unique(df1$name))) {
	df2 <- subset(df1,name==n[i])
	c <- with(df2,kruskal(value,group,p.adj = 'BH'))
	c1 <- c$groups
	c2 <- c$means
	c3 <- merge(c1[,-1,drop=F],c2,by='row.names')
	c3$name <- n[i]
	colnames(c3)[1] <- 'group'
	res <- rbind(res,c3)
	c3$group <- factor(c3$group,levels = c('Young','Old','Old+C.s'))
	df2$group <- factor(df2$group,levels = c('Young','Old','Old+C.s'))
	assign(paste0('p',i),
				 ggplot(df2,aes(x=group,y=value,fill=group))+
				 	geom_point(shape=21,size=1.8,alpha = 1,stroke=0.3,color='black',
				 						 position=position_jitterdodge(0.8),show.legend = T)+
				 	geom_bar(data=c3,stat="identity",width = 0.6,fill='transparent',color='black',lwd=0.2)+
				 	geom_errorbar(data=c3,aes(x=group,y=value,ymax = value+std,ymin = value-std),lwd=0.2,
				 								position = position_dodge(0.8),
				 								width = 0.3,color='black')+
				 	scale_fill_manual(values = c('transparent','#E24349','#3E7EB2'))+
				 	labs(y=NULL,x=NULL,title = n[i])+
				 	theme_prism(base_size = 8,	base_family = "sans",
				 							base_fontface = "plain",base_line_size = 0.2,
				 							base_rect_size = 0.3,axis_text_angle = 0,border = F )+
				 	theme(plot.title = element_text(size=8,face ='italic',color='black'),
				 				legend.position = 'none',
				 				axis.text.x = element_blank(),axis.ticks.x = element_blank(),
				 				panel.background = element_blank(),plot.background = element_blank())
	)
}
p1 <- p1+scale_y_continuous(expand = c(0,0.0001),limits = c(0,12),breaks = seq(0,12,4))
p1
p2 <- p2+scale_y_continuous(expand = c(0,0.0001),limits = c(0,3),breaks = seq(0,3,1))
p2
p3 <- p3+scale_y_continuous(expand = c(0,0.0001),limits = c(0,3),breaks = seq(0,3,1))
p3
p4 <- p4+scale_y_continuous(expand = c(0,0.0001),limits = c(0,6),breaks = seq(0,6,1.5))
p4
p5 <- p5+scale_y_continuous(expand = c(0,0.0001),limits = c(0,3),breaks = seq(0,3,1))
p5
p6 <- p6+scale_y_continuous(expand = c(0,0.0001),limits = c(0,3),breaks = seq(0,3,1))
p6
p7 <- p7+scale_y_continuous(expand = c(0,0.0001),limits = c(0,3),breaks = seq(0,3,1))
p7
p8 <- p8+scale_y_continuous(expand = c(0,0.0001),limits = c(0,3),breaks = seq(0,3,1))
p8
p9 <- p9+scale_y_continuous(expand = c(0,0.0001),limits = c(0,3),breaks = seq(0,3,1))
p9
p10 <- p10+scale_y_continuous(expand = c(0,0.0001),limits = c(0,9),breaks = seq(0,9,3))
p10
p0 <- ggarrange(p1,p2,p6,p3,p4,p5,p7,p8,p9,p10,nrow = 2,ncol=5,align = 'hv',
								common.legend = T,legend = 'bottom')
p0 <- annotate_figure(p0,left = text_grob("Relative mRNA expression", rot = 90, size = 9))
p0
#ggsave('~/Desktop/pCS aged qPCR.pdf',p0,width = 5.2,height = 2.50,dpi = 600,units = 'in')


#fig. 5D------
library(igraph);library(tidygraph);library(ggraph);library(reshape2)
ann <- read.csv("~/Desktop/Figure 5/ann_taxonmy.csv",check.names = F)
group1 <- read.csv("~/Desktop/Figure 5/groupCS.csv")
group1$group <- factor(group1$group,levels = c('Old','Old+C. s'))
path <- c("~/Desktop/Figure 5/CS 16S full network matrix/")
file_names <- c("Old.csv","Old+C. s.csv")
par(mfrow=c(1,2),mar=c(0,0,0,0),font.main=1,cex.main = 0.8,bg = 'white')
for(i in 1:length(file_names)){
	name <- gsub(".csv","",file_names[i])
	adj <- assign(name,read.csv(paste0(path,file_names[i]),row.names = 1,header = TRUE, check.names = F))
	group2 <- subset(group1,group==name)
	set.seed(123)
	g <- graph_from_adjacency_matrix(as.matrix(adj), mode = 'undirected', weight = T)
	g <- simplify(g,remove.multiple = TRUE, remove.loops = TRUE)
	g1 <- g
	node <- data.frame(degree=degree(g1),
										 cc=centr_clo(g1, mode = "all")$res,
										 bc=	centr_betw(g1, directed = FALSE)$res)
	node$bc<- (node$bc-min(node$bc))/(max(node$bc)-min(node$bc))
	node$OTU_ID <- rownames(node)
	####
	edge <- data.frame(as_edgelist(g1))  #将 igraph 默认的邻接列表格式的网络数据转为边列表格式
	edge_list <- data.frame(
		source = edge[[1]],
		target = edge[[2]])
	####
	node_list <- data.frame(
		OTU_ID = names(V(g1)),
		label = names(V(g1))
	)
	head(ann)
	node_list1 <- merge(node_list,ann,by = 'OTU_ID')
	node_list1 <- merge(node_list1,node,by='OTU_ID')
	colnames(node_list1)[1] <- 'Id' 
	node_list1$label1 <- " "
	g2 <- tbl_graph(nodes = node_list1, edges = edge_list,directed = F)
	g2 <- delete_vertices(g2, names(degree(g2)[degree(g2) == 0]))
	g2 <- simplify(g2,remove.multiple = TRUE, remove.loops = TRUE)
	V(g2)$size <- degree(g2)*0.5
	V(g2)$frame.color <- "black"  # 设置统一的边框颜色
	circle_nodes <- which(V(g2)$shape == "circle")
	non_circle_nodes <- which(V(g2)$shape != "circle")
	set.seed(123)
	layout1=layout_with_fr(g2, niter=999,grid= 'nogrid')
	plot(
		g2,
		layout = layout1,  
		vertex.label = V(g2)$label1,          
		vertex.size = V(g2)$size,   
		vertex.shape = V(g2)$shape,  
		vertex.color = "lightgray",  
		vertex.frame.color = 'black',
		vertex.frame.width = 0.2,
		edge.width = 0.2,           
		edge.arrow.size = 0.3,       
		rescale = T,
		margin = 0                  
	)
	title(main = paste0(name,'\nNodes = ',nrow(node_list),', ','edges = ',nrow(edge_list)),line = -1.8)
}
#5.5*2.5

#fig. 5E ----
df <- read.csv('~/Desktop/Figure 5/CS network attribution.csv',row.names = 1,check.names = F)
df1 <- melt(df)
df1$group <- factor(df1$group,levels = c('Old','Old+C. s'))
pnb <- ggplot(df1,aes(group,value,fill=group))+
	geom_bar(stat = 'identity',width =0.5)+
	facet_wrap(~variable,scale='free_y',nrow=3,ncol=2)+
	scale_fill_manual(values = c("Old+C. s"="#1BA3D1","Old"="#E24349"))+
	labs(x=NULL,y=NULL,title = 'Network attrubition')+
	theme_prism(base_size = 8,	base_family = "sans",
							base_fontface = "plain",base_line_size = 0.2,
							base_rect_size = 0.3,axis_text_angle = 0,border = T)+
	theme(plot.title = element_text(size=8,hjust = 0.5,vjust = 0.5,margin = margin(b=5),color = 'black'),
				strip.background = element_blank(),strip.text = element_text(size=8,hjust = 0.5,vjust = 0.5),
				panel.background = element_blank(),plot.background = element_blank(),
				legend.position = 'none')
pnb	
#ggsave('~/Desktop/pfig5e.pdf',pnb,width = 2.6,height = 4.6,dpi = 600,units = 'in')

#fig. 5F----
##keystone taxa----
library(dplyr);library(ggpubr)
deseq1 <- read.csv('~/Desktop/Figure 5/Old_node_list.csv')
deseq1 <- subset(deseq1,!cc=='NA')
deseq1$group <- 'old'
deseq2 <- read.csv('~/Desktop/Figure 5/Old+C. s_node_list.csv')
deseq2 <- subset(deseq2,!cc=='NA')
deseq2$group <- 'old+cs'
nrow(filter(deseq1,degree>=10,deseq1$cc>0.15,deseq1$bc < 0.15))
nrow(filter(deseq2,degree>=10,deseq2$cc>0.15,deseq2$bc < 0.15))
##Old----
pkey1 <- ggplot(deseq1,aes(x=cc,y=bc,size=degree))+
	geom_jitter(shape=21,fill= 'lightgray',color='black',
							alpha=0.9,stroke=0.15)+
	scale_size('Degree',range =c(0,5),breaks = c(1,5,15,25,30))+
	scale_x_continuous(expand = c(0,0.01),limits = c(0,1),breaks = seq(0,1,0.25))+
	geom_vline(xintercept = 0.15, color = 'black', size = 0.2, linetype="dashed") +
	geom_hline(yintercept = 0.15, color = 'black', size = 0.2, linetype="dashed")+
	guides(fill=guide_legend(order = 1),size=guide_legend(order = 2),shape=guide_legend(order = 3))+
	labs(x='Closeness centrality',y='Betweeness centrality',title = 'Keystone species in Old = 29')+
	theme_prism(base_size = 8,	base_family = "sans",
							base_fontface = "plain",base_line_size = 0.3,
							base_rect_size = 0.3,axis_text_angle = 0,border = T )+
	theme(plot.title = element_text(size=8,hjust = 0.5,color = 'black',margin = margin(b=5)),
				panel.background = element_blank(),plot.background = element_blank(),
				legend.position = 'bottom')
pkey1
##Old+C. S----
pkey2 <- ggplot(deseq2,aes(x=cc,y=bc,size=degree))+
	geom_point(shape=21,fill= 'lightgray',color='black',
						 alpha=0.9,stroke=0.15)+
	scale_size('Degree',range =c(0,5),breaks = c(1,5,15,25,30))+
	geom_vline(xintercept = 0.15, color = 'black', size = 0.2, linetype="dashed") +
	geom_hline(yintercept = 0.15, color = 'black', size = 0.2, linetype="dashed")+
	guides(fill=guide_legend(order = 1),size=guide_legend(order = 2),shape=guide_legend(order = 3))+
	labs(x='Closeness centrality',y='Betweeness centrality',title = 'Keystone species in Old+C. s = 39')+
	theme_prism(base_size = 8,	base_family = "sans",
							base_fontface = "plain",base_line_size = 0.2,
							base_rect_size = 0.3,axis_text_angle = 0,border = T )+
	theme(plot.title = element_text(size=8,hjust = 0.5,color = 'black',margin = margin(b=5)),
				panel.background = element_blank(),plot.background = element_blank(),
				legend.position = 'bottom')
pkey2

pkey <- ggarrange(pkey1,pkey2,nrow = 1,ncol=2,legend = 'bottom',common.legend = T)
pkey
#ggsave('~/Desktop/pfig.5f.pdf',pkey,width = 4.55,height = 2.6,dpi = 600,units = 'in')


#fig. 5G------
library(reshape2);library(vegan)
otu <- read.csv("~/Desktop/Figure 5/OTU_Old+Cs.csv",row.names = 1,check.names = F)
group <- read.csv("~/Desktop/Figure 5/groupCS.csv")
otu <- otu[which(rowSums(otu)>0),]
otu1 <- otu[,colnames(otu) %in% group$SampleID]
bray_dis <- vegdist(t(otu1), method = 'bray')      #结果以 dist 数据类型存储
nmds_dis <- metaMDS(bray_dis, k = 2)
nmds_dis$stress
nmds_dis_site <- data.frame(nmds_dis$points)
nmds_dis_site$SampleID <- rownames(nmds_dis_site)
nmds_dis_site <- merge(nmds_dis_site,group,by = 'SampleID')
nmds_dis_site <- nmds_dis_site[nmds_dis_site$SampleID %in% group$SampleID,]
nmds_dis_site$group <- factor(nmds_dis_site$group,levels = c('Old','Old+C. s'))
pNMDS <- ggplot(data = nmds_dis_site, aes(MDS1, MDS2)) +
	geom_point(aes(color = group,fill=group),width = 0.1,
							shape=21,size=2.0,stroke=0.2,alpha=1) +
	scale_color_manual(values = c('black','black'))+
	scale_fill_manual(values = c('#E24349','#3E7EB2'))+
	geom_vline(xintercept = 0, color = 'gray', size = 0.5, linetype="dashed") +
	geom_hline(yintercept = 0, color = 'gray', size = 0.5, linetype="dashed") +
	labs(x = 'NMDS1', y = 'NMDS2', title = 'NMDS plot') +
	annotate('text', label = paste('Stress =', round(nmds_dis$stress, 3)), 
					 x = -0.3, y = 0.78, size = 2.8, colour = 'black',fontface="plain")+  #标注应力函数值.
	theme_prism(base_size = 8,	base_family = "sans",
							base_fontface = "plain",base_line_size = 0.2,
							base_rect_size = 0.3,axis_text_angle = 0,border = T)+
	theme(panel.grid.major = element_blank(),
				panel.background = element_rect(color = 'black', fill = 'transparent'),
				plot.title = element_text(color='black',size=8,hjust = 0.5,margin = margin(b=5)), 
				legend.position = c(0.75,0.8),
				legend.background = element_rect(fill = "transparent",colour = "black",linewidth = 0.1),
				legend.margin = margin(0, 0, 0, 0))
pNMDS
#ggsave('~/Desktop/pfig5g.pdf',pNMDS,width = 2.2,height = 2.2,dpi = 600,units = 'in')


#fig. 5H------
plotdata <- read.csv('~/Desktop/Figure 5/deseq2 16S pac meta.csv')
plotdata$sig <- factor(plotdata$sig,levels = c('Enriched','Unchanged','Depleted'))
p1 <- ggplot(plotdata, aes(x =log2FoldChange, y =-log(padj,10))) +
	geom_vline(xintercept = c(-1, 1), color = 'black', linewidth = 0.3,linetype=2)+
	geom_hline(yintercept = -log(0.05, 10), color = 'black', linewidth= 0.3,linetype=2)+
	geom_point(aes(shape=sig,fill=sig),size=2.5,stroke=0.1,alpha=0.8) +
	scale_shape_manual(values = c(24,21,24))+
	scale_fill_manual(values  = c("#377EB8",'gray',"#E41A1C")) +
	scale_x_continuous(expand = c(0,1.5),limits = c(-30,30),breaks = seq(-30,30,15))+  
	scale_y_continuous(expand = c(0,1),limits = c(0,20),breaks = seq(-0,20,10))+  
	labs(x=bquote(Log[2]~(Foldchange)),y=bquote(-Log[10]~(italic(p)~value)),
			 title ='Differential speices analysis')+ 
	theme_prism(base_size = 8,	base_family = "sans",
							base_fontface = "plain",base_line_size = 0.3,
							base_rect_size = 0.3,axis_text_angle = 0,border = T)+
	annotate("text",label='Clostridium scindens',x=1.61957568,y=-log10(0.0004174134),size=2.6)+
	theme(panel.grid.major = element_blank(),
				panel.background = element_rect(color = 'black', fill = 'transparent'),
				plot.title = element_text(color='black',size=8,hjust = 0.5,margin = margin(b=5)), 
				legend.position = 'none')
p1
nrow(subset(plotdata,sig=='Enriched'));nrow(subset(plotdata,sig=='Depleted'))
#ggsave('~/Desktop/pfig5h.pdf',p1,width = 2.2,height = 2.2,dpi = 600,units = 'in')


#fig5I------
##CS abundance-----
te <- read.csv('~/Desktop/Figure 5/CS mean.csv');df <- read.csv('~/Desktop/Figure 5/CS abu.csv')
te$group <- factor(te$group,levels = c('Old','Old+C. s'))
df$group <- factor(df$group,levels = c('Old','Old+C. s'))
p3 <- ggplot(te, aes(x=group,y=mean)) + 
	geom_errorbar(aes(ymax = mean+sd,ymin = mean-sd),lwd=0.2,
								position = position_dodge(0.6),
								width = 0.3,color='black')+
	geom_bar(stat = 'identity',alpha=1,position = position_dodge(0.6),
					 fill='transparent',color='black',
					 width = 0.5,lwd=0.2)+
	geom_jitter(data=df,aes(x=group,y=value,fill=group),
							shape=21,size=2,alpha = 1,stroke=0.2,#color='#575F78',
							position=position_jitter(0.2),show.legend = F)+
	scale_fill_manual(values = c('#E24349','#3E7EB2'))+
	scale_y_continuous(expand = c(0,0),limits = c(0,1),breaks = seq(0,1,0.5))+
	labs(x=NULL,y='Relative abundance (%)',title = 'C. scindens')+
	theme_prism(base_size = 8,	base_family = "sans",
							base_fontface = "plain",base_line_size = 0.2,
							base_rect_size = 0.3,axis_text_angle = 0,border = F)+
	theme(panel.grid.major = element_blank(),
				strip.text = element_text(size=8,hjust = 0.5),
				panel.background = element_blank(),plot.background = element_blank(),
				plot.title = element_text(color='black',face = 'italic',size=8,hjust = 0.5,margin = margin(b=5)), 
				legend.position = 'none')+
	geom_signif(data=df,aes(x=group,y=value),
							comparisons = list(c('Old','Old+C. s')),
							map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
							y_position = 0.8,
							tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
							test = "t.test",test.args = 'two.sided')
p3
#ggsave('~/Desktop/pfig5i_CS abu.pdf',p3,width = 1.3,height = 1.6,dpi = 600,units = 'in')

## CS IAA-----
library(ggh4x)
df <- read.csv('~/Desktop/Figure 5/CS colon and feces IAA change.csv')
df1 <- melt(df,id.vars = c("SampleID","group","tissue"),variable.name = 'name',value.name = 'value')
te1 = aggregate(df1$value,                      
								by  =list(df1$group,df1$name,df1$tissue),FUN='mean')
te2 = aggregate(df1$value,                      
								by  =list(df1$group,df1$name,df1$tissue),FUN='sd')
te <- cbind(te1,te2[4])
colnames(te) <- c('group','name',"tissue",'mean','sd')
te$group <- factor(te$group,levels = c("Young",'Old','Old+C.s'))
df1$group <- factor(df1$group,levels = c("Young",'Old','Old+C.s'))
piaa1 <- ggplot(te, aes(x=group,y=mean)) + 
	geom_errorbar(aes(ymax = mean+sd,ymin = mean-sd),lwd=0.2,
								position = position_dodge(0.6),
								width = 0.3,color='black')+
	geom_bar(stat = 'identity',alpha=1,position = position_dodge(0.6),
					 fill='transparent',color='black',
					 width = 0.5,lwd=0.2)+
	geom_jitter(data=df1,aes(x=group,y=value,fill=group),
							shape=21,size=2,alpha = 1,stroke=0.2,#color='#575F78',
							position=position_jitter(0.2),show.legend = F)+
	scale_fill_manual(values = c('transparent','#E24349','#3E7EB2'))+
	facet_wrap(~tissue, nrow=1,scale='free_y')+
	facetted_pos_scales(
		y = list(tissue == 'Colon' ~ 
						 	scale_y_continuous(expand = c(0,0), limits = c(0,8),breaks=seq(0,8,2)),
						 tissue == 'Feces' ~ 
						 	scale_y_continuous(expand = c(0,0), limits = c(0,80),breaks=seq(0,80,20))))+
	labs(x=NULL,y='Concentration (ng/mg)',title =NULL)+
	theme_prism(base_size = 8,	base_family = "sans",
							base_fontface = "plain",base_line_size = 0.2,
							base_rect_size = 0.3,axis_text_angle = 0,border = F)+
	theme(plot.title = element_text(size=8,hjust = 0.5,vjust = 0.5,margin = margin(b=5),color = 'black'),
				strip.background = element_blank(),strip.text = element_text(size=8,hjust = 0.5,vjust = 0.5),
				panel.background = element_blank(),plot.background = element_blank(),
				legend.position = 'none')+
	geom_signif(data=df1,aes(x=group,y=value),
							comparisons = list(c("Young",'Old'),c('Old',"Old+C.s")),
							map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),y_position = 5,step_increase = 0.1,
							tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
							test = "t.test",test.args = 'two.sided')
piaa1
#
#ggsave('~/Desktop/pfig5i_IAA.pdf',piaa1,width = 2.5,height = 1.50,dpi = 600,units = 'in')

#Finally, these plots are arranged into fig. 5 in AI software and further modify legends.
