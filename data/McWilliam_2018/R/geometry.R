
###################################### 
# geometric equations

hemisphere<-function(d, tb, ...){
SA<-(2*pi*(d/2)^2)
TB<-SA*tb
VOL<-((2/3)*pi*(d/2)^3)
data.frame(SA, TB, VOL)}

plate2d<-function(d, th, tb, ...){
SA<-pi*((d/2)^2)
TB<-SA*tb
VOL<-SA*th
data.frame(SA, TB, VOL)}

up_plates<-function(d, br, bpa, th, tb, ...){
SA<-(pi*(d/2)^2*(bpa*( 2*pi*((br)^2)) ))
TB<-SA*tb
VOL<-th*(SA)
data.frame(SA,TB,VOL)}

open_cone<-function(d, bh, th, tb, ...){
SA<-(2*pi*(d/2)*sqrt((d/2)+bh))
TB<-SA*tb
VOL<-th*(SA/2)
data.frame(SA,TB, VOL)}

branches<-function(d, bh, br, bpa, tb, ...){ 
SA<-(pi*((d/2)^2)*(bpa*(2*pi*br*bh)+(pi*(br^2))))
TB<-SA*tb
VOL<-(pi*((d/2)^2)*(bpa*(pi*(br^2)*bh)))
data.frame(SA,TB, VOL)}

plate2d_proj<-function(d, bh, br, bpa, th, tb, ...){
SA<-(pi*((d/2)^2))+(pi*((d/2)^2)*(bpa*((2*br*bh)+(pi*(br^2)))))
TB<-SA*tb
VOL<-(pi*((d/2)^2)*th)+(pi*((d/2)^2)*(bpa*(pi*(br^2)*bh)))
data.frame(SA,TB, VOL)}
