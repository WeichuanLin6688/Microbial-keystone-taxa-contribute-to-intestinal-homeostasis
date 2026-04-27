#fig. 6 ----
#fig. 6B ----
##IAA AB-PAS count-------
library(ggplot2);library(agricolae);library(ggprism)
df <- read.csv("~/Desktop/Figure 6/IAA AB-PAS count.csv",check.names = F)
df$group <- factor(df$group,levels = c('Young','Old','Old+IAA'))
c <- with(df,kruskal(value,group,p.adj = 'BH'))
c1 <- c$groups
c2 <- c$means
c3 <- merge(c1[,-1,drop=F],c2,by='row.names')
colnames(c3)[1] <- 'group'
c3$group <- factor(c3$group,levels = c('Young','Old','Old+IAA'))
pabpas <- ggplot(df,aes(x=group,y=value,fill=group))+
	geom_point(shape=21,size=1.8,alpha = 1,stroke=0.3,color='black',
						 position=position_jitterdodge(0.8),show.legend = F)+
	geom_bar(data=c3,stat="identity",width = 0.60,fill='transparent',color='black',lwd=0.2)+
	geom_errorbar(data=c3,aes(x=group,y=value,ymax = value+std,ymin = value-std),lwd=0.2,
								position = position_dodge(0.8),
								width = 0.3,color='black')+
	scale_fill_manual(values = c('transparent','#E24349','#FC7F28'))+
	labs(y="Cell counts\nper millimeter crypt",x=NULL,title = 'Goblet cells')+
	theme_prism(base_size = 8,	base_family = "sans",
							base_fontface = "plain",base_line_size = 0.2,
							base_rect_size = 0.3,axis_text_angle = 0,border = F )+
	theme(plot.title = element_text(size=9,face ='plain',color='black',margin = margin(b=5)),
				legend.position = 'right',
				axis.text.x = element_text(size=8,angle = 45,hjust = 1,vjust = 1),
				panel.background = element_blank(),plot.background = element_blank())		+
	geom_signif(comparisons = list(c('Young','Old'),c('Old','Old+IAA')),
							map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
							y_position = c(65,75),
							tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
							test = "t.test",test.args = 'two.sided')+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,90),breaks = seq(0,90,30))
pabpas 
#ggsave('~/Desktop/pfig6b_IAA AB-PAS count.pdf',pabpas,width = 1.2,height = 1.5,dpi = 600,units = 'in')

##IAA colon γH2AX----
df <- read.csv("~/Desktop/Figure 6/IAA γH2AX.csv",check.names = F)
df$group <- factor(df$group,levels = c('Young','Old','Old+IAA'))
c <- with(df,kruskal(mean,group,p.adj = 'BH'))
c1 <- c$groups
c2 <- c$means
c3 <- merge(c1[,-1,drop=F],c2,by='row.names')
colnames(c3)[1] <- 'group'
c3$group <- factor(c3$group,levels = c('Young','Old','Old+IAA'))
df$group <- factor(df$group,levels = c('Young','Old','Old+IAA'))
with(df[df$group %in% c('Old','Old+IAA') ,],t.test(mean~group))
with(df[df$group %in% c('Young','Old') ,],t.test(mean~group))
pH2AX <- ggplot(df,aes(x=group,y=mean,fill=group))+
	geom_point(shape=21,size=2,alpha = 1,stroke=0.2,color='black',
						 position=position_jitterdodge(0.8),show.legend = T)+
	geom_bar(data=c3,stat="identity",width = 0.6,fill='transparent',color='black',lwd=0.2)+
	geom_errorbar(data=c3,aes(x=group,y=mean,ymax = mean+std,ymin = mean-std),lwd=0.2,
								position = position_dodge(0.8),
								width = 0.3,color='black')+
	scale_fill_manual(values = c('transparent','#E24349','#FC7F28'))+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,3),breaks = seq(0,3,1))+
	labs(y="Area (%)",x=NULL,title = 'γH2A.X')+
	theme_prism(base_size = 8,	base_family = "sans",
							base_fontface = "plain",base_line_size = 0.2,
							base_rect_size = 0.3,axis_text_angle = 0,border = F )+
	theme(plot.title = element_text(size=8,face ='plain',color='black',margin = margin(b=5)),
				legend.position = 'none',
				axis.text.x = element_text(size=8,angle = 45,hjust = 1,vjust = 1),
				panel.background = element_blank(),plot.background = element_blank())+
	geom_signif(data=df,aes(x=group,y=mean),
							comparisons = list(c('Young','Old'),c('Old','Old+IAA')),
							map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
							y_position = c(1.9,2.05),
							tip_length = 0,color='black',textsize = 3,vjust = -0.2,size = 0.2,
							test = "t.test",test.args = 'two.sided')
pH2AX
#ggsave('~/Desktop/pfig6b_IAA pH2AX.pdf',pH2AX,width = 1.0,height = 1.4,dpi = 600,units = 'in')


##IAA colon CLDN10----
df <- read.csv("~/Desktop/Figure 6/IAA CLDN10.csv",check.names = F)
df$group <- factor(df$group,levels = c('Young','Old','Old+IAA'))
c <- with(df,kruskal(mean,group,p.adj = 'BH'))
c1 <- c$groups
c2 <- c$means
c3 <- merge(c1[,-1,drop=F],c2,by='row.names')
colnames(c3)[1] <- 'group'
c3$group <- factor(c3$group,levels = c('Young','Old','Old+IAA'))
df$group <- factor(df$group,levels = c('Young','Old','Old+IAA'))
with(df[df$group %in% c('Old','Old+IAA') ,],t.test(mean~group))
with(df[df$group %in% c('Young','Old') ,],t.test(mean~group))
pcldn10 <- ggplot(df,aes(x=group,y=mean,fill=group))+
	geom_point(shape=21,size=2,alpha = 1,stroke=0.2,color='black',
						 position=position_jitterdodge(0.8),show.legend = T)+
	geom_bar(data=c3,stat="identity",width = 0.6,fill='transparent',color='black',lwd=0.2)+
	geom_errorbar(data=c3,aes(x=group,y=mean,ymax = mean+std,ymin = mean-std),lwd=0.2,
								position = position_dodge(0.8),
								width = 0.3,color='black')+
	scale_fill_manual(values = c('transparent','#E24349','#FC7F28'))+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,9),breaks = seq(0,9,3))+
	labs(y="Area (%)",x=NULL,title = 'CLDN10')+
	theme_prism(base_size = 8,	base_family = "sans",
							base_fontface = "plain",base_line_size = 0.2,
							base_rect_size = 0.3,axis_text_angle = 0,border = F )+
	theme(plot.title = element_text(size=8,face ='plain',color='black',margin = margin(b=5)),
				legend.position = 'none',
				axis.text.x = element_text(size=8,angle = 45,hjust = 1,vjust = 1),
				#axis.text.x = element_blank(),axis.ticks.x = element_blank(),
				panel.background = element_blank(),plot.background = element_blank())+
	geom_signif(data=df,aes(x=group,y=mean),
							comparisons = list(c('Young','Old'),c('Old','Old+IAA')),
							map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
							y_position = c(7.5,8.2),
							tip_length = 0,color='black',textsize = 3,vjust = -0.2,size = 0.2,
							test = "t.test",test.args = 'two.sided')
pcldn10
#ggsave('~/Desktop/pfig6b_IAA cldn10.pdf',pcldn10,width = 1.0,height = 1.4,dpi = 600,units = 'in')


#fig. 6C IAA 肠道PCR------
df1 <- read.csv("~/Desktop/Figure 6/IAA colon PCR.csv",check.names = F)
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
	c3$group <- factor(c3$group,levels = c('Young','Old','Old+IAA'))
	df2$group <- factor(df2$group,levels = c('Young','Old','Old+IAA'))
	assign(paste0('p',i),
				 ggplot(df2,aes(x=group,y=value,fill=group))+
				 	geom_point(shape=21,size=1.8,alpha = 1,stroke=0.2,color='black',
				 						 position=position_jitterdodge(0.8),show.legend = T)+
				 	geom_bar(data=c3,stat="identity",width = 0.65,fill='transparent',color='black',lwd=0.2)+
				 	geom_errorbar(data=c3,aes(x=group,y=value,ymax = value+std,ymin = value-std),lwd=0.2,
				 								position = position_dodge(0.8),
				 								width = 0.35,color='black')+
				 	scale_fill_manual(values = c('transparent','#E24349','#FC7F28'))+
				 	labs(y=NULL,x=NULL,title = n[i])+
				 	theme_prism(base_size = 8,	base_family = "sans",
				 							base_fontface = "plain",base_line_size = 0.2,
				 							base_rect_size = 0.3,axis_text_angle = 0,border = F )+
				 	theme(plot.title = element_text(size=8,face ='italic',color='black',margin = margin(b=5)),
				 				legend.position = 'top',
				 				axis.text.x = element_text(size=7,angle = 20,hjust = 0.9,vjust = 0.9),
				 				panel.background = element_blank(),plot.background = element_blank())
	)
}
p1 <- p1+	geom_signif(comparisons = list(c('Young','Old'),c('Old','Old+IAA')),
											map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
											y_position = c(7,7.5),
											tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
											test = "t.test",test.args = 'two.sided')+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,9),breaks = seq(0,9,3))

p2 <- p2+	geom_signif(comparisons = list(c('Young','Old'),c('Old','Old+IAA')),
											map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
											y_position = c(4,5),
											tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
											test = "t.test",test.args = 'two.sided')+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,9),breaks = seq(0,9,3))

p3 <- p3+	geom_signif(comparisons = list(c('Young','Old'),c('Old','Old+IAA')),
											map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
											y_position = c(3,4),
											tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
											test = "t.test",test.args = 'two.sided')+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,6),breaks = seq(0,6,2))

p4 <- p4+	geom_signif(comparisons = list(c('Young','Old'),c('Old','Old+IAA')),
											map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
											y_position =  c(3,4),
											tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
											test = "t.test",test.args = 'two.sided')+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,6),breaks = seq(0,6,2))

p5 <- p5+	geom_signif(comparisons = list(c('Young','Old'),c('Old','Old+IAA')),
											map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
											y_position = c(3,4),
											tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
											test = "t.test",test.args = 'two.sided')+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,6),breaks = seq(0,6,2))

p6 <- p6 +geom_signif(comparisons = list(c('Young','Old'),c('Old','Old+IAA')),
											map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
											y_position =  c(3,4),
											tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
											test = "t.test",test.args = 'two.sided')+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,6),breaks = seq(0,6,2))

p7 <- p7 +geom_signif(comparisons = list(c('Young','Old'),c('Old','Old+IAA')),
											map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
											y_position = c(3,5),
											tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
											test = "t.test",test.args = 'two.sided')+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,6),breaks = seq(0,6,2))

p8 <- p8 +geom_signif(comparisons = list(c('Young','Old'),c('Old','Old+IAA')),
											map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
											y_position = c(3,4),
											tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
											test = "t.test",test.args = 'two.sided')+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,6),breaks = seq(0,6,2))

p0 <- ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,nrow = 2,ncol=4,align = 'hv',
								common.legend = T,legend = 'none'
)
p0 <- annotate_figure(p0,left = text_grob("Relative mRNA expression", rot = 90, size = 8))
p0
#ggsave('~/Desktop/pfig6c_IAA PCR.pdf',p0,width = 3.8,height = 2.7,dpi = 600,units = 'in')

#fig. 6E-----
##IAA AhR AB-PAS count------
df2 <- read.csv("~/Desktop/Figure 6/IAA AhR AB-PAS.csv",check.names = F)
c <- with(df2,kruskal(value,group,p.adj = 'BH'))
c1 <- c$groups
c2 <- c$means
c3 <- merge(c1[,-1,drop=F],c2,by='row.names')
colnames(c3)[1] <- 'group'
c3$group <- factor(c3$group,levels = c('OldAhrfl/fl','OldAhrfl/fl+IAA','OldAhrIEC','OldAhrIEC+IAA'))
df2$group <- factor(df2$group,levels = c('OldAhrfl/fl','OldAhrfl/fl+IAA','OldAhrIEC','OldAhrIEC+IAA'))
p2 <- ggplot(df2,aes(x=group,y=value,fill=group))+
	geom_point(shape=21,size=1.8,alpha = 1,stroke=0.2,color='black',
						 position=position_jitterdodge(0.8),show.legend = T)+
	geom_signif(data=df2,aes(x=group,y=value),
							comparisons = list(c('OldAhrfl/fl','OldAhrfl/fl+IAA'),
																 c('OldAhrfl/fl','OldAhrIEC'),
																 c('OldAhrIEC','OldAhrIEC+IAA')),
							map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
							y_position = c(70,80,75),
							tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
							test = "t.test",test.args = 'two.sided')+
	geom_bar(data=c3,stat="identity",width = 0.6,fill='transparent',color='black',lwd=0.2)+
	geom_errorbar(data=c3,aes(x=group,y=value,ymax = value+std,ymin = value-std),lwd=0.2,
								position = position_dodge(0.8),
								width = 0.3,color='black')+
	scale_fill_manual(values = c('#E24349','#FC7F28',"#999999",'#3B7EB5'))+
	labs(y="Cell counts\nper millimeter crypt",x=NULL,title = 'Goblet cells')+
	theme_prism(base_size = 8,	base_family = "sans",
							base_fontface = "plain",base_line_size = 0.2,
							base_rect_size = 0.3,axis_text_angle = 0,border = F )+
	theme(plot.title = element_text(size=8,face ='plain',color='black',margin = margin(b=5)),
				legend.position = 'none',
				axis.text.x = element_text(size=7,angle = 45,hjust = 1,vjust = 1),
				panel.background = element_blank(),plot.background = element_blank())+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,90),breaks = seq(0,90,30))
p2
#ggsave('~/Desktop/pfig6e_AhR AB-PAS.pdf',p2,width = 1.4,height = 1.7,dpi = 600,units = 'in')

##IAA Ahr γH2AX------
df2 <- read.csv("~/Desktop/Figure 6/IAA AhR γH2AX.csv",check.names = F)
df2$value <- df2$value
c <- with(df2,kruskal(value,group,p.adj = 'BH'))
c1 <- c$groups
c2 <- c$means
c3 <- merge(c1[,-1,drop=F],c2,by='row.names')
colnames(c3)[1] <- 'group'
c3$group <- factor(c3$group,levels = c('OldAhrfl/fl','OldAhrfl/fl+IAA','OldAhrIEC','OldAhrIEC+IAA'))
df2$group <- factor(df2$group,levels = c('OldAhrfl/fl','OldAhrfl/fl+IAA','OldAhrIEC','OldAhrIEC+IAA'))
p1 <- ggplot(df2,aes(x=group,y=value,fill=group))+
	geom_point(shape=21,size=1.8,alpha = 1,stroke=0.2,color='black',
						 position=position_jitterdodge(0.8),show.legend = T)+
	geom_signif(data=df2,aes(x=group,y=value),
							comparisons = list(c('OldAhrfl/fl','OldAhrfl/fl+IAA'),
																 c('OldAhrfl/fl','OldAhrIEC'),
																 c('OldAhrIEC','OldAhrIEC+IAA')),
							map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
							y_position = c(1.5,2,2.5),
							tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
							test = "t.test",test.args = 'two.sided')+
	geom_bar(data=c3,stat="identity",width = 0.6,fill='transparent',color='black',lwd=0.2)+
	geom_errorbar(data=c3,aes(x=group,y=value,ymax = value+std,ymin = value-std),lwd=0.2,
								position = position_dodge(0.8),
								width = 0.3,color='black')+
	scale_fill_manual(values = c('#E24349','#FC7F28',"#999999",'#3B7EB5'))+
	labs(y="Area (%)",x=NULL,title = 'γH2A.X')+
	theme_prism(base_size = 8,	base_family = "sans",
							base_fontface = "plain",base_line_size = 0.2,
							base_rect_size = 0.3,axis_text_angle = 0,border = F )+
	theme(plot.title = element_text(size=8,face ='plain',color='black',margin = margin(b=5)),
				legend.position = 'none',
				axis.text.x = element_text(size=7,angle = 45,hjust = 1,vjust = 1),
				panel.background = element_blank(),plot.background = element_blank())+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,3),breaks = seq(0,3,1))
p1
#ggsave('~/Desktop/pfig6e_IAA AhR γh2ax.pdf',p1,width = 1.2,height = 1.7,dpi = 600,units = 'in')

#fig. 6F-----
df1 <- read.csv("~/Desktop/Figure 6/IAA AhR PCR.csv")
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
	c3$group <- factor(c3$group,levels = c('OldAhrfl/fl','OldAhrfl/fl+IAA','OldAhrIEC','OldAhrIEC+IAA'))
	df2$group <- factor(df2$group,levels = c('OldAhrfl/fl','OldAhrfl/fl+IAA','OldAhrIEC','OldAhrIEC+IAA'))
	assign(paste0('p',i),
				 ggplot(df2,aes(x=group,y=value,fill=group))+
				 	geom_point(shape=21,size=1.8,alpha = 1,stroke=0.2,color='black',
				 						 position=position_jitterdodge(0.8),show.legend = T)+
				 	geom_bar(data=c3,stat="identity",width = 0.65,fill='transparent',color='black',lwd=0.2)+
				 	geom_errorbar(data=c3,aes(x=group,y=value,ymax = value+std,ymin = value-std),lwd=0.2,
				 								position = position_dodge(0.8),
				 								width = 0.35,color='black')+
				 	scale_fill_manual(values = c('#E24349','#FC7F28',"#999999",'#3B7EB5'))+
				 	labs(y=NULL,x=NULL,title = n[i])+
				 	theme_prism(base_size = 8,	base_family = "sans",
				 							base_fontface = "plain",base_line_size = 0.2,
				 							base_rect_size = 0.3,axis_text_angle = 0,border = F )+
				 	theme(plot.title = element_text(size=8,face ='italic',color='black',margin = margin(b=5)),
				 				legend.position = 'rigjht',
				 				axis.text.x = element_text(size=7,angle = 45,hjust = 1,vjust = 1),
				 				panel.background = element_blank(),plot.background = element_blank())
	)
}
p1 <- p1+	geom_signif(comparisons = list(c('OldAhrfl/fl','OldAhrfl/fl+IAA'),
																				 c('OldAhrfl/fl','OldAhrIEC'),
																				 c('OldAhrIEC','OldAhrIEC+IAA')),
											map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
											y_position = c(2,6,7),
											tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
											test = "t.test",test.args = 'two.sided')+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,9),breaks = seq(0,9,3))

p2 <- p2+	geom_signif(comparisons = list(c('OldAhrfl/fl','OldAhrfl/fl+IAA'),
																				 c('OldAhrfl/fl','OldAhrIEC'),
																				 c('OldAhrIEC','OldAhrIEC+IAA')),
											map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
											y_position = c(1.5,2,2.5),
											tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
											test = "t.test",test.args = 'two.sided')+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,9),breaks = seq(0,9,3))

p3 <- p3+	geom_signif(comparisons = list(c('OldAhrfl/fl','OldAhrfl/fl+IAA'),
																				 c('OldAhrfl/fl','OldAhrIEC'),
																				 c('OldAhrIEC','OldAhrIEC+IAA')),
											map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
											y_position = c(2.5,3.5,4.5),
											tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
											test = "t.test",test.args = 'two.sided')+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,6),breaks = seq(0,6,2))

p4 <- p4+	geom_signif(comparisons = list(c('OldAhrfl/fl','OldAhrfl/fl+IAA'),
																				 c('OldAhrfl/fl','OldAhrIEC'),
																				 c('OldAhrIEC','OldAhrIEC+IAA')),
											map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
											y_position = c(2.5,3.5,4.5),
											tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
											test = "t.test",test.args = 'two.sided')+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,6),breaks = seq(0,6,2))

p5 <- p5+	geom_signif(comparisons = list(c('OldAhrfl/fl','OldAhrfl/fl+IAA'),
																				 c('OldAhrfl/fl','OldAhrIEC'),
																				 c('OldAhrIEC','OldAhrIEC+IAA')),
											map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
											y_position = c(2.5,3.5,4.5),
											tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
											test = "t.test",test.args = 'two.sided')+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,6),breaks = seq(0,6,2))

p6 <- p6 +geom_signif(comparisons = list(c('OldAhrfl/fl','OldAhrfl/fl+IAA'),
																				 c('OldAhrfl/fl','OldAhrIEC'),
																				 c('OldAhrIEC','OldAhrIEC+IAA')),
											map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
											y_position = c(2.5,3.5,4.5),
											tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
											test = "t.test",test.args = 'two.sided')+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,6),breaks = seq(0,6,2))

p7 <- p7 + geom_signif(comparisons = list(c('OldAhrfl/fl','OldAhrfl/fl+IAA'),
																					c('OldAhrfl/fl','OldAhrIEC'),
																					c('OldAhrIEC','OldAhrIEC+IAA')),
											 map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
											 y_position = c(4.5,5.2,2.5),
											 tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
											 test = "t.test",test.args = 'two.sided')+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,6),breaks = seq(0,6,2))

p8 <- p8 + geom_signif(comparisons = list(c('OldAhrfl/fl','OldAhrfl/fl+IAA'),
																					c('OldAhrfl/fl','OldAhrIEC'),
																					c('OldAhrIEC','OldAhrIEC+IAA')),
											 map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
											 y_position = c(2.5,3.5,1.5),
											 tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
											 test = "t.test",test.args = 'two.sided')+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,6),breaks = seq(0,6,2))

p0 <- ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,nrow = 2,ncol=4,align = 'hv',
								common.legend = T,legend = 'none'
)
p0 <- annotate_figure(p0,left = text_grob("Relative mRNA expression", rot = 90, size = 8))
p0
#ggsave('~/Desktop/pfig6f_IAA Ahr PCR.pdf',p0,width = 4,height = 3.5,dpi = 600,units = 'in')


#fig. 6G IAA  transcriptome heatmap-----
df1 <- read.csv('~/Desktop/Figure 6/AhR transcriptome.csv')
n <- unique(df1$name)
res <- NULL
for (i in 1: length(n)) {
	df2 <- subset(df1,name==n[i])
	df2$value1 <- scale(df2$value)
	c <- with(df2,kruskal(value,group,group = T,p.adj = 'BH'))
	df2$pvalue <- c$statistics$p.chisq
	res <- rbind(res,df2)
}
res$SampleID <- factor(res$SampleID,levels = rev(c('A1','A2','A3',
																									 'B1','B2','B3',
																									 'C1','C2','C3')))
res$name <- factor(res$name,levels = c(unique(df1$name)))
p1 <- ggplot(res,aes(name,SampleID))+
	geom_tile(aes(fill=value1),lwd=0.1,
						color='lightgray',width=0.95, height=0.95,alpha=1)+
	scale_fill_gradient2('Gene counts\n(z-score)',limits=c(min(res$value1), max(res$value1)),
											 low='#3B7EB5',mid = 'transparent',high='#E41A1C')+
	guides(fill=guide_colorbar(barwidth=0.4,barheight=4))+
	labs(y='Ahrfl/fl        Ahrfl/fl+IAA        AhrIEC')+
	theme(plot.background = element_blank(),panel.background = element_blank(),
				plot.title = element_text(size=8,hjust = 0.5,vjust = 0.5),
				axis.title.y = element_text(size=8,hjust = 0.5,vjust = 0.5),
				axis.text.y  = element_blank(),
				axis.ticks.x = element_line(size=0.3),
				axis.ticks.y = element_blank(),legend.position = 'right',
				#axis.ticks.y = element_line(size=0.3),
				axis.text.x = element_text(size=8,angle = 30,vjust = 1,hjust = 1,color='black'),
				#axis.text.y = element_text(size=8,vjust = 0.3,hjust = 1,color='black'),
				legend.title = element_text(size=8,hjust = 0,vjust = 0.5),
				axis.title.x = element_blank())+
	annotate("segment", x = 0.2, xend = 0.2, y = 0.7, yend = 3.3,size = 1,color='black')+
	annotate("segment", x = 0.2, xend = 0.2, y = 3.7, yend = 6.3,size = 1,color='black')+
	annotate("segment", x = 0.2, xend = 0.2, y = 6.7, yend = 9.3,size = 1,color='black')
p1
#ggsave('~/Desktop/pfig6g_AhR transcriptome.pdf',p1,width = 4.55,height = 2.35,dpi = 600,units = 'in')

#fig. 6H-----
res <- read.csv('~/Desktop/Figure 6/IAA AhR transcriptome PCR.csv')
te1 = aggregate(res$value1,                      
								by  = list(res$name,res$group),FUN='mean')
te2 = aggregate(res$value1,                      
								by  = list(res$name,res$group),FUN='sd')
te <- cbind(te1,te2[3])
names(te) = c('name','group','mean','sd')
te$name <- factor(te$name,levels = unique(res$name))
te$group <- factor(te$group,levels = c('Young','Old','Old+IAA'))
res$group <- factor(res$group,levels = c('Young','Old','Old+IAA'))
p <- ggplot(te, aes(x=name,y=mean,fill=group,group)) + 
	geom_hline(yintercept = 1, color = 'black', size = 0.3, linetype="dashed") +
	geom_errorbar(aes(ymax = mean+sd,ymin = mean-sd),lwd=0.2,
								position = position_dodge(0.8),
								width = 0.3,color='black')+
	geom_bar(stat = 'identity',alpha=0.1,position = position_dodge(0.8),color='black',
					 show.legend = F,
					 width = 0.6,lwd=0.2)+
	scale_fill_manual(values = c('transparent','#E24349','#FC7F28'))+
	geom_jitter(data=res,aes(x=name,y=value1,color=group,group = group),
							shape=21,size=1.2,alpha = 0.8,stroke=0.2,color='black',#fill='gray',
							position=position_jitterdodge(0.2),show.legend = T)+
	scale_color_manual(values = c('transparent','#E24349','#FC7F28'))+
	scale_y_continuous(expand = c(0,0),limits =  c(0,6.0),breaks = seq(0,6,2))+
	labs(x=NULL,y='Relative mRNA expression',title=NULL)+
	guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))+
	theme_prism(base_size = 8,	base_family = "sans",
							base_fontface = "plain",base_line_size = 0.2,
							base_rect_size = 0.3,axis_text_angle = 0,border = F)+
	theme(legend.position = c(0.5,0.95),legend.direction = 'horizontal',
				plot.title = element_text(size=8,hjust = 0.5,margin = margin(b=5)),
				legend.text = element_text(size=8,hjust = 0.5,vjust = 0.5),
				axis.text.x = element_text(face = 'italic',size=8,hjust = 1,vjust = 1,angle = 30))
p
#ggsave('~/Desktop/pfig6h.pdf',p,width = 4.55,height = 2.05,dpi = 600,units = 'in')

#Finally, these plots are arranged into fig. 6 in AI software and further modify legends.


