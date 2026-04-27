#fig. 4 ----
#fig. 4B ----
library(ggplot2);library(gggenes)
genes1 <- data.frame(
	gene = c('','','','',''),
	color = c('','','aldh','',''),
	start = c(541,1553,3001,4564,5840),  
	end = c(1494,2902,4482,5838,6482), 
	strand = c("+","+","+","+","+")   
)
pa1 <- ggplot(genes1, aes(xmin = start, xmax = end, y = gene, label = color,
													fill = color, forward = strand == "+")) +
	geom_gene_arrow(size=0.2,
									arrowhead_height = unit(6, "mm"), arrowhead_width = unit(1, "mm")
	)+
	geom_gene_label(size=10) +
	scale_x_continuous(breaks = c(3001,4482),labels = c(2225797,"2227278 bp"))+
	labs(title = expression(italic('C. scindens')~' genome GCA_004295125.1'),x=NULL,y=NULL)+
	scale_fill_manual(values = c('gray','#EB8678','#8BB0D0','#90D3C7'))+
	theme_genes()+theme(panel.grid.major.y = element_line(colour = "grey", 
																												linewidth = 1.5),
											plot.title = element_text(size=9,color='black',hjust = 0.5),
											panel.background = element_blank(),
											legend.position = 'none',plot.background = element_blank(),
											axis.line.x = element_line(colour = "black", linewidth = 0.3), 
											axis.ticks.x = element_line(colour = "black", size=0.3),
											axis.text.x = element_text(size=9,color='black'))
print(pa1)

genes2 <- data.frame(
	gene = c('','','','',''),
	color = c('','','amie','',''),
	start = c(2299290,2299642,2300384,2302005,2302912
	), 
	end = c(2299484,2300247,2301004,2302568,2304246
	),   
	strand = c("+","+","+","+","+"
	)  
)
pa2 <- ggplot(genes2, aes(xmin = start, xmax = end, y = gene, label = color,
													fill = color, forward = strand == "+")) +
	geom_gene_arrow(size=0.2,
									arrowhead_height = unit(6, "mm"), arrowhead_width = unit(1, "mm")
	)+
	geom_gene_label(size=10) +
	scale_x_continuous(breaks = c(2300384,2301004),labels = c(2300384,"2301004 bp"))+
	labs(title = NULL,x=NULL,y=NULL)+
	scale_fill_manual(values = c('gray','#8BB0D0'))+
	theme_genes()+theme(panel.grid.major.y = element_line(colour = "grey", 
																												linewidth = 1.5),
											plot.title = element_text(size=9,color='black',hjust = 0.5),
											panel.background = element_blank(),
											legend.position = 'none',plot.background = element_blank(),
											axis.line.x = element_line(colour = "black", linewidth = 0.3), 
											axis.ticks.x = element_line(colour = "black", size=0.3),
											axis.text.x = element_text(size=9,color='black'))
print(pa2)

p8 <- ggarrange(pa1,pa2,align = 'hv',ncol = 1,nrow=2,heights = c(1,0.9))
p8
#ggsave('~/Desktop/pfig4b.pdf',p8,width = 3,height = 2.2,dpi = 600,units = 'in')

#fig. 4D ----
library(ggplot2);library(ggprism);library(splines);library(ggtrendline)
data <- read.csv('~/Desktop/Figure 4/CS growth rate.csv')
x <- as.numeric(data$time)
y <- as.numeric(data$value)
p1 <- ggtrendline(x,y,model = "line3P",
									linecolor = '#FC7F28',linetype = 1,linewidth = 0.6,
									CI.level = 0.95,CI.fill ='#E5E5E5',CI.alpha = 1,
									CI.color = "black",CI.lty = 0,CI.lwd = 0.3,
									summary = TRUE,show.eq = TRUE,yhat = FALSE,eq.x = 40,eq.y = 0.48,
									Rname = 1,Pname = 0,rrp.x = 30,rrp.y = 0.43,
									text.col = "black",eDigit = 3,eSize = 3.0,xlab = NULL,ylab = NULL)+
	labs(x='Time (hours)',y=expression("OD"[600~'nm']),
			 title = expression('The growth rate of C. scindens'))+
	geom_jitter(aes(x,y),shape=21,size=2.5,#position = position_dodge(0.8),
							height = 0.01,width = 0,
							stroke=0.3,color='black')+
	scale_y_continuous(expand = c(0,0),limits = c(0.2,0.5),breaks = seq(0.2,0.5,0.1))+
	scale_x_continuous(limits = c(0,72),breaks = c(0,12,24,48,72))+
	theme_prism(base_size = 9,	base_family = "sans",
							base_fontface = "plain",base_line_size = 0.3,
							base_rect_size = 0.3,axis_text_angle = 0,border = FALSE)+
	theme(legend.position = 'none',
				plot.background = element_blank(),panel.background = element_blank(),
				plot.title = element_text(size=9,hjust = 0.5,margin = margin(b=5)))
p1
#ggsave('~/Desktop/pfig4d.pdf',p1,width = 2.5,height = 2.3,units = 'in',dpi = 600)

#fig. 4F ----
df <- read.csv("~/Desktop/Figure 4/CS IAA.csv",check.names = F)
df1 <- melt(df,id.vars = 'group',variable.name = 'name',value.name = 'value')
n <- unique(df1$name)
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
	c3$group <- factor(c3$group,levels = c("CK","C. s"))
	df2$group <- factor(df2$group,levels =c("CK","C. s"))
	assign(paste0('p',i),
				 ggplot(df2,aes(x=group,y=value,color=group))+
				 	geom_bar(data=c3,aes(x=group,y=value,fill=group),stat="identity",alpha=0.9,
				 					 width = 0.5,color='black',lwd=0.2)+
				 	geom_errorbar(data=c3,aes(x=group,y=value,ymax = value+std,ymin = value-std),lwd=0.2,
				 								position = position_dodge(0.5),
				 								width = 0.25,color='black')+
				 	scale_fill_manual(values = c('lightgray','#FC7F28'))+
				 	geom_point(shape=21,size=2.,alpha = 1,stroke=0.3,color='black',fill='lightgray',
				 						 position=position_jitter(0.2),show.legend = F)+
				 	labs(y=NULL,x=NULL,title = n[i])+
				 	theme_prism(base_size = 8,	base_family = "sans",
				 							base_fontface = "plain",base_line_size = 0.2,
				 							base_rect_size = 0.3,axis_text_angle = 0,border = FALSE )+
				 	theme(plot.title = element_text(size=9,face ='plain',color='black',margin = margin(b=5)),
				 				legend.position = 'none',
				 				panel.background = element_blank(),plot.background = element_blank())
	)
}
p1 <- p1+geom_signif(comparisons = list(c("CK","C. s")),
										 map_signif_level=function(p) sprintf("p = %.3f", p),
										 step_increase = 0.6,y_position = 650,
										 tip_length = 0,color='black',textsize = 3.0,vjust = -0.8,size = 0.3,
										 test = "t.test",test.args = 'two.sided')+
	scale_y_continuous(expand = c(0,0),limits = c(0,800),breaks = seq(0,800,200))
p2 <- p2+geom_signif(comparisons = list(c("CK","C. s")),
										 map_signif_level=function(p) sprintf("p = %.3f", p),
										 step_increase = 0.6,y_position = 32000,
										 tip_length = 0,color='black',textsize = 3.0,vjust = -0.8,size = 0.3,
										 test = "t.test",test.args = 'two.sided')+
	scale_y_continuous(expand = c(0,0),limits = c(0,40000),breaks = seq(0,40000,10000))
p3 <- p3+geom_signif(comparisons = list(c("CK","C. s")),
										 map_signif_level=function(p) sprintf("p = %.3f", p),
										 step_increase = 0.6,y_position = 32000,
										 tip_length = 0,color='black',textsize = 3.0,vjust = -0.8,size = 0.3,
										 test = "t.test",test.args = 'two.sided')+
	scale_y_continuous(expand = c(0,0),limits = c(0,40000),breaks = seq(0,40000,10000))

p8 <- ggarrange(p1,p2,p3,align = 'hv',ncol = 3,nrow=1)
p9 <- annotate_figure(p8,left = text_grob("Concentration (nmol/mL)", rot = 90, size = 9))
p9
#ggsave('~/Desktop/pfig4f.pdf',p9,width = 4.5,height = 1.8,dpi = 600,units = 'in')

#Finally, these plots are arranged into fig. 4 in AI software and further modify legends.




