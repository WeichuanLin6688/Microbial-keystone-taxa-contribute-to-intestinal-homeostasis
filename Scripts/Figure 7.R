#fig. 7 ----
#fig. 7B ----
##IAA AhR敲除 细胞PCR------
df1 <- read.csv("~/Desktop/Figure 7/Cell AhR IAA.csv",check.names = F)
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
	c3$group <- factor(c3$group,levels = c('DOXO','DOXO+IAA','shAhr+DOXO','shAhr+DOXO+IAA'))
	df2$group <- factor(df2$group,levels = c('DOXO','DOXO+IAA','shAhr+DOXO','shAhr+DOXO+IAA'))
	assign(paste0('p',i),
				 ggplot(df2,aes(x=group,y=value,fill=group))+
				 	geom_point(shape=21,size=1.8,alpha = 1,stroke=0.2,color='black',
				 						 position=position_jitterdodge(0.8),show.legend = T)+
				 	geom_bar(data=c3,stat="identity",width = 0.65,fill='transparent',color='black',lwd=0.2)+
				 	geom_errorbar(data=c3,aes(x=group,y=value,ymax = value+std,ymin = value-std),lwd=0.2,
				 								position = position_dodge(0.8),
				 								width = 0.35,color='black')+
				 	scale_fill_manual(values = c('#E24349','#FC7F28',"#999999",'#3B7EB5'))+
				 	#scale_y_continuous(expand = c(0,0.0001),limits = c(0,9),breaks = seq(0,9,3))+
				 	labs(y=NULL,x=NULL,title = n[i])+
				 	theme_prism(base_size = 8,	base_family = "sans",
				 							base_fontface = "plain",base_line_size = 0.2,
				 							base_rect_size = 0.3,axis_text_angle = 0,border = F )+
				 	theme(plot.title = element_text(size=8,face ='italic',color='black',margin = margin(b=5)),
				 				legend.position = 'rigjht',
				 				#axis.text.x = element_blank(),axis.ticks.x = element_blank(),
				 				#axis.text.x = element_text(size=8,angle = 0,hjust = 0.5,vjust = 0.5),
				 				axis.text.x = element_text(size=7,angle = 30,hjust = 0.9,vjust = 0.9),
				 				panel.background = element_blank(),plot.background = element_blank())
	)
}
p1 <- p1+	geom_signif(comparisons = list(c('DOXO','DOXO+IAA'),
																				 c('DOXO','shAhr+DOXO'),
																				 c('shAhr+DOXO','shAhr+DOXO+IAA')),
											map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
											y_position = c(3,6,4.5),
											tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
											test = "t.test",test.args = 'two.sided')+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,9),breaks = seq(0,9,3))

p2 <- p2+	geom_signif(comparisons = list(c('DOXO','DOXO+IAA'),
																				 c('DOXO','shAhr+DOXO'),
																				 c('shAhr+DOXO','shAhr+DOXO+IAA')),
											map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
											y_position = c(3,5,3.5),
											tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
											test = "t.test",test.args = 'two.sided')+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,9),breaks = seq(0,9,3))

p3 <- p3+	geom_signif(comparisons = list(c('DOXO','DOXO+IAA'),
																				 c('DOXO','shAhr+DOXO'),
																				 c('shAhr+DOXO','shAhr+DOXO+IAA')),
											map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
											y_position = c(2.5,3.5,2),
											tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
											test = "t.test",test.args = 'two.sided')+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,6),breaks = seq(0,6,2))

p4 <- p4+	geom_signif(comparisons = list(c('DOXO','DOXO+IAA'),
																				 c('DOXO','shAhr+DOXO'),
																				 c('shAhr+DOXO','shAhr+DOXO+IAA')),
											map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
											y_position = c(2.5,3.5,2),
											tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
											test = "t.test",test.args = 'two.sided')+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,6),breaks = seq(0,6,2))

p5 <- p5+	geom_signif(comparisons = list(c('DOXO','DOXO+IAA'),
																				 c('DOXO','shAhr+DOXO'),
																				 c('shAhr+DOXO','shAhr+DOXO+IAA')),
											map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
											y_position = c(4.5,5.2,3),
											tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
											test = "t.test",test.args = 'two.sided')+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,6),breaks = seq(0,6,2))

p6 <- p6 + geom_signif(comparisons = list(c('DOXO','DOXO+IAA'),
																					c('DOXO','shAhr+DOXO'),
																					c('shAhr+DOXO','shAhr+DOXO+IAA')),
											 map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
											 y_position = c(3,4,3),
											 tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
											 test = "t.test",test.args = 'two.sided')+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,6),breaks = seq(0,6,2))

p7 <- p7 + geom_signif(comparisons = list(c('DOXO','DOXO+IAA'),
																					c('DOXO','shAhr+DOXO'),
																					c('shAhr+DOXO','shAhr+DOXO+IAA')),
											 map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
											 y_position = c(3,4.5,3),
											 tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
											 test = "t.test",test.args = 'two.sided')+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,6),breaks = seq(0,6,2))

p8 <- p8 + geom_signif(comparisons = list(c('DOXO','DOXO+IAA'),
																					c('DOXO','shAhr+DOXO'),
																					c('shAhr+DOXO','shAhr+DOXO+IAA')),
											 map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
											 y_position = c(2.5,4,2.5),
											 tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
											 test = "t.test",test.args = 'two.sided')+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,6),breaks = seq(0,6,2))

p0 <- ggarrange(p1,p2,p8,p6,p4,p7,p5,p3,nrow = 2,ncol=4,align = 'hv',
								common.legend = T,legend = 'none')
p0 <- annotate_figure(p0,left = text_grob("Relative mRNA expression", rot = 90, size = 8))
p0
#ggsave('~/Desktop/pfig6b.pdf',p0,width = 8,height = 1.65,dpi = 600,units = 'in')

#fig7. C-----
##Cell IAA AhR γH2AX------
df2 <- read.csv("~/Desktop/Figure 7/Cell AhR γH2AX.csv",check.names = F)
c <- with(df2,kruskal(value,group,p.adj = 'BH'))
c1 <- c$groups
c2 <- c$means
c3 <- merge(c1[,-1,drop=F],c2,by='row.names')
colnames(c3)[1] <- 'group'
c3$group <- factor(c3$group,levels = c('DOXO','DOXO+IAA','shAhr+DOXO','shAhr+DOXO+IAA'))
df2$group <- factor(df2$group,levels = c('DOXO','DOXO+IAA','shAhr+DOXO','shAhr+DOXO+IAA'))
p1 <- ggplot(df2,aes(x=group,y=value,fill=group))+
	geom_point(shape=21,size=1.8,alpha = 1,stroke=0.2,color='black',
						 position=position_jitterdodge(0.8),show.legend = T)+
	geom_signif(data=df2,aes(x=group,y=value),
							comparisons = list(c('DOXO','DOXO+IAA'),
																 c('DOXO','shAhr+DOXO'),
																 c('shAhr+DOXO','shAhr+DOXO+IAA')),
							map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
							y_position = c(6,7,6.5),
							tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
							test = "t.test",test.args = 'two.sided')+
	geom_bar(data=c3,stat="identity",width = 0.65,fill='transparent',color='black',lwd=0.2)+
	geom_errorbar(data=c3,aes(x=group,y=value,ymax = value+std,ymin = value-std),lwd=0.2,
								position = position_dodge(0.8),
								width = 0.35,color='black')+
	scale_fill_manual(values = c('#E24349','#FC7F28',"#999999",'#3B7EB5'))+
	labs(y="Percentage (%)",x=NULL,title = 'γH2AX')+
	theme_prism(base_size = 8,	base_family = "sans",
							base_fontface = "plain",base_line_size = 0.2,
							base_rect_size = 0.3,axis_text_angle = 0,border = F )+
	theme(plot.title = element_text(size=8,face ='plain',color='black',margin = margin(b=5)),
				legend.position = 'none',
				axis.text.x = element_text(size=7,angle = 30,hjust = 0.9,vjust = 0.9),
				panel.background = element_blank(),plot.background = element_blank())+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,9),breaks = seq(0,9,3))
p1
#ggsave('~/Desktop/pfig6c_Cell AhR γh2ax.pdf',p1,width = 1.2,height = 1.6,dpi = 600,units = 'in')

##Cell IAA AhR CLDN10------
df2 <- read.csv("~/Desktop/Figure 7/Cell AhR CLDN10.csv",check.names = F)
c <- with(df2,kruskal(value,group,p.adj = 'BH'))
c1 <- c$groups
c2 <- c$means
c3 <- merge(c1[,-1,drop=F],c2,by='row.names')
colnames(c3)[1] <- 'group'
c3$group <- factor(c3$group,levels = c('DOXO','DOXO+IAA','shAhr+DOXO','shAhr+DOXO+IAA'))
df2$group <- factor(df2$group,levels = c('DOXO','DOXO+IAA','shAhr+DOXO','shAhr+DOXO+IAA'))
p2 <- ggplot(df2,aes(x=group,y=value,fill=group))+
	geom_point(shape=21,size=1.8,alpha = 1,stroke=0.2,color='black',
						 position=position_jitterdodge(0.8),show.legend = T)+
	geom_signif(data=df2,aes(x=group,y=value),
							comparisons = list(c('DOXO','DOXO+IAA'),
																 c('DOXO','shAhr+DOXO'),
																 c('shAhr+DOXO','shAhr+DOXO+IAA')),
							map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
							y_position = c(1.3,1.5,0.3),
							tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
							test = "t.test",test.args = 'two.sided')+
	geom_bar(data=c3,stat="identity",width = 0.65,fill='transparent',color='black',lwd=0.2)+
	geom_errorbar(data=c3,aes(x=group,y=value,ymax = value+std,ymin = value-std),lwd=0.2,
								position = position_dodge(0.8),
								width = 0.35,color='black')+
	scale_fill_manual(values = c('#E24349','#FC7F28',"#999999",'#3B7EB5'))+
	labs(y="Percentage (%)",x=NULL,title = 'CLDN10')+
	theme_prism(base_size = 8,	base_family = "sans",
							base_fontface = "plain",base_line_size = 0.2,
							base_rect_size = 0.3,axis_text_angle = 0,border = F )+
	theme(plot.title = element_text(size=8,face ='plain',color='black',margin = margin(b=5)),
				legend.position = 'none',
				axis.text.x = element_text(size=7,angle = 30,hjust = 0.9,vjust = 0.9),
				panel.background = element_blank(),plot.background = element_blank())+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,2),breaks = seq(0,2,1))
p2
#ggsave('~/Desktop/pfig6_Cell AhR CLDN10.pdf',p2,width = 1.2,height = 1.6,dpi = 600,units = 'in')

#fig. 7F--------
df3 <- read.csv("~/Desktop/Figure 7/IAA AhR Cldn10 promoter.csv",check.names = F)
group <- df3[!duplicated(df3$group),]
c <- with(df3,kruskal(value,group,p.adj = 'BH'))
c1 <- c$groups
c2 <- c$means
c3 <- merge(c1[,-1,drop=F],c2,by='row.names')
colnames(c3)[1] <- 'group'
c3 <- merge(c3,group[,-2])
c3$group <- factor(c3$group,levels = unique(df3$group))
df3$group <- factor(df3$group,levels = unique(df3$group))
p3 <- ggplot(data=c3,aes(x=group,y=value,fill=treatment))+
	geom_bar(data=c3,aes(x=group,y=value,fill=treatment),stat="identity",#fill="transparent",
					 width = 0.5,color='black',lwd=0.2,alpha=1)+
	geom_errorbar(data=c3,aes(x=group,y=value,ymax = value+std,ymin = value-std),lwd=0.2,
								position = position_dodge(0.8),
								width = 0.25,color='black')+
	scale_fill_manual(values = c('transparent',"transparent",'transparent'))+
	geom_point(data=df3,aes(x=group,y=value,color=treatment,fill=treatment),
						 shape=1,size=1.8,alpha = 1,stroke=0.3,fill='transparent',
						 position=position_jitterdodge(0.8),show.legend = T)+
	scale_color_manual(values = c('#FC7F28','#3B7EB5',"#999999"))+
	labs(y="Relative luciferase activity\n(FL/RL)",x=NULL,title = NULL)+
	theme_prism(base_size = 8,	base_family = "sans",
							base_fontface = "plain",base_line_size = 0.2,
							base_rect_size = 0.3,axis_text_angle = 0,border = F )+
	theme(plot.title = element_text(size=8,face ='plain',color='black',margin = margin(b=5)),
				legend.position = 'none',
				axis.text.x = element_text(size=7,angle = 30,hjust = 0.9,vjust = 0.9),
				panel.background = element_blank(),plot.background = element_blank())+
	scale_y_continuous(expand = c(0,0.0001),limits = c(0,15),breaks = seq(0,15,5))+
	geom_signif(data=df3,aes(x=group,y=value),
							comparisons = list(c('Cldn10 promotor2','Cldn10 promotor3'),
																 c('Cldn10 promotor3','Cldn10XRE promotor5'),
																 c('Cldn10 promotor3','Cldn10 promotor7')
							),
							map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),
							y_position = c(11,13,11.5),
							tip_length = 0,color='black',textsize = 2.8,vjust = -0.2,size = 0.2,
							test = "t.test",test.args = 'two.sided')
p3
#ggsave('~/Desktop/pfig6f_IAA AhR cldn10 promoter.pdf',p3,width = 3.0,height = 2.6,dpi = 600,units = 'in')

#Finally, these plots are arranged into fig. 7 in AI software and further modify legends.
