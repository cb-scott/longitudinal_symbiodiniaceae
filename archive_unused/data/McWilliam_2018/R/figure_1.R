
# vectors

theme1 <- theme_void()+
theme(panel.border=element_rect(size=0.5, fill=NA, colour="black"), 
strip.text=element_text(size=7), 
plot.margin = unit(c(1,1,1,1), "mm"), 
plot.title=element_text(size=7, hjust=.5))

global <- ggplot(data=NULL, aes(PC1, PC2))+
  geom_polygon(data=axes[chull(axes$PC1, axes$PC2),], col="grey", fill=NA) + 
  geom_point(data=axes, col="grey",size=0.2)+theme1

vec1 <- global+ 
  geom_segment(data=vecs,aes(0,0,xend=PC1*5.2,yend=-PC2*5.2), col="red", size=0.3)+
  geom_text(data=vecs, aes(PC1*6, -PC2*5.5, label=label), size=2.5,fontface="bold")+
  lims(x=c(-4.8, 5.16), y=c(-3.6, 3.6))+
  labs(x=paste(vars[1], "%"), y=paste(vars[2], "%"), title="Vectors (PC1-2)")+ 
  theme(axis.title.y=element_text(size=7, angle=90), axis.title.x=element_text(size=7))

vec2 <- ggplot(data=NULL, aes(PC3, PC4))+
  geom_polygon(data=axes[chull(axes$PC3, axes$PC4),], col="grey", fill=NA) +
  geom_point(data=axes, colour="grey", size=0.2)+
  geom_segment(data=vecs,aes(0,0,xend=PC3*5.2,yend=-PC4*5.2), col="red", size=0.3)+
  geom_text(data=vecs, aes(PC3*6, -PC4*5.5, label=label), size=2.5,fontface="bold")+
  lims(x=c(-4.8, 5.16), y=c(-4.8, 3.6))+
  labs(x=paste(vars[3], "%"), y=paste(vars[4], "%"), title="Vectors (PC3-4)")+  theme1 + theme(axis.title.y=element_text(size=7, angle=90), axis.title.x=element_text(size=7))

# convex hulls

# order
ord <- info[order(info$order),"region"]
hulls$region<-factor(hulls$region, levels=ord)
points$region<-factor(points$region, levels=ord) 

# group
points$g <- info$g[match(points$region, info$region)]
hulls$g <- info$g[match(hulls$region, info$region)]

# label
vols <- unique(hulls[,c("region","n","vol", "g")])

pp <- function(x, col){
  global+
  geom_point(data=axes, colour="grey", size=.8)+
    geom_polygon(data=hulls[hulls$g==x,], aes(fill=region), col="black", alpha = 0.6)+
    geom_point(data=points[points$g==x,], aes(fill=region), size=0.8, shape=21, stroke=0.25)+
    geom_text(data=vols[vols$g==x,], aes(x=-2.7, y=-2.9, label=paste(round(vol,1), "%")), size=2.5, col="slategrey")+
    geom_text(data=vols[vols$g==x,], size=2.5, col="slategrey",aes(x=2.8, y=-2.9, label=paste("S =",round(n,1))))+ facet_wrap(~region, nrow=1)+
    lims(x=c(-4, 4.3), y=c(-3,3))+
    scale_fill_manual(values=cols)+guides(fill="none")+
    theme(plot.background=element_rect(colour=col, size=1))
  }
  
# plot

fig1 <- plot_grid(pp("a", "slategrey"), map,
          plot_grid(plot_grid(vec1, vec2),pp("b", "slategrey"), pp("c", "red"), nrow=1, rel_widths=c(0.66,0.66,1), scale=0.98), 
          nrow=3, rel_heights=c(0.9,1,0.9), scale=0.99)
          
          
          
          

# FIGURE S1

hulls34 <- NULL
for (r in info$region){
df <- axes[!is.na(pres_abs[,r]),]
hulls34 <-rbind(hulls34, cbind(df[chull(df$PC3, df$PC4),], region=r))
}

hulls34$region <- factor(hulls34$region, levels=ord)
hulls34$g <- info$g[match(hulls34$region, info$region)]

pp2 <- function(x, col){
ggplot()+
geom_point(data=axes, aes(PC3, PC4), col="grey", size=0.8)+
geom_polygon(data=axes[chull(axes$PC3, axes$PC4),], aes(PC3, PC4), col="grey", fill=NA)+
geom_polygon(data=hulls34[hulls34$g==x,], aes(PC3, PC4, fill=region), col="black", alpha = 0.6)+
geom_point(data=points[points$g==x,], aes(PC3, PC4, fill=region), size=0.8, shape=21, stroke=0.25)+
facet_wrap(~region, nrow=1)+
geom_text(data=vols[vols$g==x,], aes(x=-2.7, y=-2.9, label=paste(round(vol,1), "%")), size=2.5, col="slategrey")+
geom_text(data=vols[vols$g==x,], size=2.5, col="slategrey",aes(x=2.8, y=-2.9, label=paste("S =",round(n,1))))+ facet_wrap(~region, nrow=1)+
    lims(x=c(-4, 4.3), y=c(-3,3))+
    theme1+
    scale_fill_manual(values=cols)+
    guides(fill="none")+
    theme(plot.background=element_rect(colour=col, size=1))
    }
          
fig_s1 <- plot_grid(pp2("a", "slategrey"),
          plot_grid(plot_grid(vec1, vec2),pp2("b", "slategrey"), pp2("c", "red"), nrow=1, rel_widths=c(0.66,0.66,1), scale=0.98), 
          nrow=2, rel_heights=c(0.9,1,0.9), scale=0.99)
          
     fig_s1     
          
          


