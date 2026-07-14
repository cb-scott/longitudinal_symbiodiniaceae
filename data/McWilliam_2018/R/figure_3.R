
# figure 3

theme3 <- theme_bw()+
theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
axis.text=element_blank(),  legend.position=c(0.82,0.91), axis.text.y=element_blank(), panel.grid = element_blank(), 
axis.ticks=element_blank(), legend.text=element_text(size=6), legend.title=element_blank(),legend.key.size=unit(2,"mm"),legend.key = element_blank(),legend.background=element_rect(colour="grey"),plot.title = element_text(size=10))

columnar<-readPNG("data/img/columnar.png")
columnar<-rasterGrob(columnar, interpolate=TRUE, )
digitate<-readPNG("data/img/digitate.png")
digitate<-rasterGrob(digitate, interpolate=TRUE)
encrusting<-readPNG("data/img/encrusting.png")
encrusting<-rasterGrob(encrusting, interpolate=TRUE)
hispidose<-readPNG("data/img/hispidose.png")
hispidose<-rasterGrob(hispidose, interpolate=TRUE)
laminar<-readPNG("data/img/laminar.png")
laminar<-rasterGrob(laminar, interpolate=TRUE)
massive<-readPNG("data/img/massive.png")
massive<-rasterGrob(massive, interpolate=TRUE)
branching<-readPNG("data/img/branching.png")
branching<-rasterGrob(branching, interpolate=TRUE)
solitary<-readPNG("data/img/solitary.png")
solitary<-rasterGrob(solitary, interpolate=TRUE)
submassive<-readPNG("data/img/submassive.png")
submassive<-rasterGrob(submassive, interpolate=TRUE)
tabular<-readPNG("data/img/tabular.png")
tabular<-rasterGrob(tabular, interpolate=TRUE)

axes$repro <- data[,"reproductive_mode"]
axes$repro[axes$repro=="MIXED"] <- "brooder"
doms <- axes[!is.na(axes$region),]

A <- ggplot(doms, aes(x=PC1, y=PC2))+
        stat_density2d(geom="polygon", bins=4, alpha=0.25, fill="slategrey")+
        stat_density2d(colour="black", bins=4, size=0.3)+
        geom_polygon(data = axes[chull(axes$PC1, axes$PC2),], colour="grey", fill=NA, alpha = 0.5) +
        geom_point(size=0.5, colour="slategrey")+
        geom_point(aes(fill=region, shape=region), size=2, stroke=0.25)+
        lims(x=c(-4, 4.3), y=c(-3,3))+
        scale_fill_manual(values=c("#ff7f00","turquoise"))+
        scale_shape_manual(values=c(24,21))+
        annotation_custom(columnar, xmin=1, xmax=2.5, ymax=-1.8, ymin=-3)+
        annotation_custom(branching, xmin=2.9, xmax=4.6, ymin=-1, ymax=-1.8)+
        annotation_custom(hispidose, xmin=2, xmax=3.7, ymin=-0.9, ymax=-2.6)+
        annotation_custom(digitate, xmin=1.7, xmax=3.3, ymin=0.9, ymax=2.4)+
        annotation_custom(tabular, xmin=2.5, xmax=4.3, ymin=0.4, ymax=2)+
        annotation_custom(laminar, xmin=-2.8, xmax=-0.5, ymin=2, ymax=3.6)+
        annotation_custom(encrusting, xmin=-4.3, xmax=-2.9, ymin=0.5, ymax=2.1)+
        annotation_custom(solitary, xmin=-4.4, xmax=-3, ymin=-0, ymax=-1.5)+
        annotation_custom(massive, xmin=-4, xmax=-2, ymin=-1.4, ymax=-3.1)+
        annotation_custom(submassive, xmin=-2.5, xmax=-1.3, ymin=-1.7, ymax=-3.1)+
        theme3

panelb <- list()
for (d in unique(doms$region)){
        plot <- ggplot(doms[doms$region==d,], aes(PC1,PC2))+
                stat_density2d(geom="polygon", aes(fill=..level..), alpha=0.25, bins=200)+
                scale_fill_distiller(palette="Spectral", name="Density",limits=c(0,0.105), breaks=0.005*seq(0,10,by=2))+
                guides(color="none",fill="none")+
                ggtitle(d)+theme3
        panelb[[d]] <- plot
        }

B <- plot_grid(
        panelb[["Great Barrier Reef"]]+lims(y=c(-3,3.5), x=c(-4.8,5.4))+
                stat_density2d(colour="black", alpha=1, bins=5.8),
        panelb[["Caribbean"]]+ lims(y=c(-4.2,4.2), x=c(-4.7,5.1))+
                stat_density2d(colour="black", alpha=1, bins=6.2)
        )

cc <- cbind(rep(unique(doms$region),each=2),rep(c("Biomass", "Accretion"),2))
nbin <-  10

panelc <- list()
for(i in 1:nrow(cc)){
        df <- subset(funs, region==cc[i,1])
        df$funct <- df[,cc[i,2]] 
        plot<- ggplot()+
                stat_summary_2d(data=df, aes(x=PC1, y=PC2, z=1),bins=nbin, size=0.15, colour="black", fill="grey")+
                stat_summary_2d(data=df[!is.na(df$funct),], aes(x=PC1, y=PC2, z=log(funct)),fun=mean, bins=nbin, size=0.15, colour="black")+
                guides(fill="none")+theme3+lims(y=c(-3,3.5), x=c(-4.8,5.4))+
                annotate("text",x = Inf, y = Inf, hjust = 1.1, vjust=1.6, label=cc[i,2], size=2)
        means<-ggplot_build(plot)$data[[2]]$value
        pp <- plot + scale_fill_distiller(palette="Spectral",limits=c(min(means), max(means)+1))
        panelc[[i]] <- pp
        }

C <- plot_grid(
        plot_grid(panelc[[1]], NULL, panelc[[2]], nrow=1, rel_widths=c(1,-0.09,1)),
        plot_grid(panelc[[3]]+lims(y=c(-4.2,4.2), x=c(-4.7,5.1)), NULL,
                  panelc[[4]]+lims(y=c(-4.2,4.2), x=c(-4.7,5.1)), 
                  nrow=1, rel_widths=c(1,-0.09,1)), nrow=1
        )

paneld <- list()
for (d in unique(doms$region)){
        df <- doms[doms$region==d,]
        pd <- ggplot(data=clust, aes(PC1, PC2))+ 
                stat_ellipse(aes(group=factor(k)), geom="polygon",level=0.95, alpha=0.1, colour="grey", size=0.5)+
                geom_point(data=df[!is.na(df$repro),], aes(fill=repro, shape=repro), stroke=0.3, size=1.3)+
                lims(x=c(-4, 4.3), y=c(-3,3.2))+
                scale_shape_manual(values=c(24,21), labels=c("Brooders", "Spawners"))+
                scale_fill_manual(values=c("red","grey"), labels=c("Brooders", "Spawners"))+
                theme3+
                theme(legend.margin=margin(1,1,1,1), legend.position=c(0.82,0.88), legend.text=element_text(size=5))
        paneld[[d]] <- pd
        }

D <-  plot_grid(paneld[["Great Barrier Reef"]]+guides(fill="none", shape="none"),
                paneld[["Caribbean"]])

fig3 <- plot_grid(A, NULL, B, C, D, 
                  nrow=5, rel_heights=c(1,-0.03, 0.53, 0.28, 0.46), 
                  labels=c("A","","B","C","D"), 
                  label_size=10, label_x=c(0.03,0,0.03,0.03,0.03), 
                  label_y=c(0.98,0,0.84,0.95, 0.95)
                  )
fig3