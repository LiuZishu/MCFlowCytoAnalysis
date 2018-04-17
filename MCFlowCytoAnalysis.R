#1. Introduction======================================================================
#The R-script used for the article: 
#"Neutral mechanisms and niche differentiation in steady-state insular microbial communities reveald by single cell analysis"
#Liu et al., 2018
#----------------------------------------------------------------------------------
#Required packages
library(vegan)
library(reshape2)
library(ggplot2)
library(grid)
#----------------------------------------------------------------------------------
#Data sets:
#abundance.txt: relative cells abundance per subcommunity per sample
#Explination for each colunm
#"Sample_name": name of samples;
#"Group": reactors where samples picked from;
#"Time_d": sampling time (day)
#"G1-68": subcommunities
abundance<-read.table("abundance.txt",header=T,sep="\t")

#abio.txt:abiotic parameters
#Explination for each colunm
#"Group": reactors where samples picked from;
#"Time": sampling time (day)
# Abiotic parameters, such as "Tmp": temperature; "pH":pH; "EC": electric conductivity; "NH4": ammonium concertration;...
abio=read.table("abio.txt",header=TRUE,row.names=1,na.strings = "NA",blank.lines.skip = FALSE)

#All.txt:relative cells abundance per subcommunity and corresponding abiotic parameters
All=read.table("All.txt",header=T,sep="\t")

#TestRefSpace.txt:the test result of whether reference space was built per reactor
EcoSat=read.table("TestRefSpace.txt",header=T)

#According to reactor groups, D(q=0) of dominant subcommunities were divided into two files:"diversity.txt" for controls; "diversityDs.txt" for disturbed reactors.
divCs=read.table("diversity.txt",header=T)
divDs=read.table("diversityDs.txt",header=T)

#te.txt: Temperature data for disturbed reactors (D1,D2 and D3)
te=read.table("temperature.txt",header=T)
#----------------------------------------------------------------------------------


#2. Codes==================================================================================
#2.1 Diversity calculation
#The function SC.div was built to calculate gate based Hill numbers
#In the function, requirements are q: the order(0, 1, 2) for calculation; tr: threshold for distinguish dominant subcommunities (i.e. 1/68=0.01471[1.47%] in this study); data: data set used for calculation.

SC.div<-function(q,tr,data){
  abundantdata=data.frame(data)
  nor.data=data.frame(data)
  D=matrix(nrow=nrow(data),ncol=1)
  rownames(D)=rownames(data)
  if(tr==0){
    nor.data=data
  }else
  {
    for(i in 1:ncol(data)) 
    {
      for(j in 1:nrow(data)) 
      {
        if(data[j,i]>=tr){
          abundantdata[j,i]=data[j,i] 
        }else{
          abundantdata[j,i]=0
        }
      }
    }
    
    for(j in 1:nrow(data)) 
    {
      rowsum=sum(abundantdata[j,])
      for(i in 1:ncol(data)) 
      {
        nor.data[j,i]=abundantdata[j,i]/rowsum
      }
    }
  }
  caldata=nor.data
  if (q==1){
    D[,1]=exp(diversity(caldata,index="shannon"))
  }else{
    for(i in 1:ncol(data)) 
    {
      for(j in 1:nrow(data)) 
      {
        if(nor.data[j,i]==0){
          caldata[j,i]=0
        }else{
          caldata[j,i]=nor.data[j,i]^q
        }
      }
    }
    
    for(k in 1:nrow(data))
    {
      D[k,1]=(sum(caldata[k,]))^(1/(1-q))
    }
  }
  
  return(D)  
}
#For example, to calculate the Dq=0 of dominant subcommunities in this study
D0dSc=SC.div(q=0,tr=0.01471,data=abundance[4:71])
rownames(D0dSc)=abundance$Sample_name
#Dq=0 of all subcommunities in this study
D0allSc=SC.div(q=0,tr=0,data=abundance[4:71])
#Dq=1 of all subcommunities in this study
D1allSc=SC.div(q=1,tr=0,data=abundance[4:71])
#Dq=2 of all subcommunities in this study
D2allSc=SC.div(q=2,tr=0,data=abundance[4:71])

#D(q=0,1,2) for all subcommunities and D(q=0) of dominant subcommunities were calculated and compared in S9.
AlphaDiversity=data.frame('Sample'=abundance$Sample_name,'D0dSc'=D0dSc,'D0allSc'=D0allSc,'D1allSc'=D1allSc, 'D2allSc'=D2allSc)
View(AlphaDiversity)
#According to reactor groups, D(q=0) of dominant subcommunities were divided into two files:"diversity.txt" for controls; "diversityDs.txt" for disturbed reactors.

#============================================================================
#2.2 Testing of reference space
#2.2.1 Last 10 samples per reactor used to test the reference space
#e.g. data of C1

data=abundance[2:66,]
rownames(data)=data[,1]
numberOfGates = 68#68 subcommunities in this study
data$referencePhase = FALSE
data[56:65, "referencePhase"] <- TRUE
referenceState <- colMeans(data[data$referencePhase==TRUE, 4:(numberOfGates + 3)])
data$Canberra <- apply(data[4:(numberOfGates+3)], 1, function(x) dist(rbind(referenceState, x), method = "canberra"))
data$Canberra <- data$Canberra / numberOfGates
r_Canberra <- max(data$Canberra[data$referencePhase==TRUE]) 
#The boundary of reference space is determined as r_Canberra, which is 0.23 for C1.
#Results as r_Canberra and values of "data$Canberra" per reactor write out as the file EcoSat.txt, which will be used to show the test result for the reference space per reactor (see 2.3.1 and 2.3.3).


#============================================================================
#2.3 Codes for the main text
#2.3.1 Figure 1

#Command for creating pdf file
pdf("Figure_1.pdf",height=11.3,width=10)
layout(matrix(c(1,2,3,4,5,5), ncol=2, byrow=TRUE), heights=c(5,5,1))
par(mai=c(0.55,0.55,0.25,0.25),omi=c(0.1,0.1,0.2,0.1),xpd=FALSE)
#------------------------------------------------------------------
#Virsulization of the reference space
plot(x=EcoSat$time,y=EcoSat$C1,type="n",pch=19,col="gray50",ylim=c(0,1),xaxt="n",xlab="Time (d)",ylab="Deviation from reference state",font.lab=2,cex.lab=1.2)

abline(h=0.23,col="gray50",lty=2,lwd=2)
abline(h=0.24,col="grey",lty=2,lwd=2)
points(x=EcoSat$time,y=EcoSat$C1,type="o",pch=19,col="gray50",cex=1.5,lwd=2)
points(x=EcoSat$time,y=EcoSat$C2,type="o",pch=19,col="grey",cex=1.5,lwd=2)

axis(side=1,at= seq(0,91,7))
mtext("A", side = 3, line = -1,adj=0,outer = TRUE,cex=1.5)
#------------------------------------------------------------------
#Dissimilarity analysis with relative cells abundance per subcommunities (all five reactor)
rownames(abundance)<-abundance[,1]
bcnmds = metaMDS(abundance[1:326,4:71], distance="bray", autotransform=FALSE, zerodist="add",k=2,try=500,trymax=1000)
mds.out=bcnmds

#parameters for visulization
size=1+3.5*abundance$Time_d[2:66]/91
#points size increasing with sampling time
group=c("Inoculum",rep("Control",130),rep("Disturbed",195))
reactor=as.character(abundance$Group)
time=abundance$Time_d[2:66]
col=c("gray50","grey","deepskyblue2","brown3","darkseagreen3")


plot(mds.out, type="n",font.lab=2,cex.lab=1.2,bg="transparent")
mtext("B", side = 3, line = -1,adj=0.5,outer = TRUE,cex=1.5)
# mark reference state
arrows(mds.out$points[1,1],mds.out$points[1,2],mds.out$points[2,1],mds.out$points[2,2], length = 0.08, angle = 15, lwd=1.5, col ="gray50")
arrows(mds.out$points[1,1],mds.out$points[1,2],mds.out$points[67,1],mds.out$points[67,2], length = 0.08, angle = 15, lwd=1.5, col ="grey")
points(mds.out$points[1,1],mds.out$points[1,2], col ="black", pch = 21, cex = 1.5,bg="black",lwd=3)
points(mds.out$points[67:131,], col ="grey", type="l",lty=2)
points(mds.out$points[67:131,], col ="grey", pch = 21,bg="white",cex = size,type="p",lwd=2)
points(mds.out$points[2:66,], col ="gray50", type="l",lty=2)
points(mds.out$points[2:66,], col ="gray50", pch = 21,bg="white",cex = size,type="p",lwd=2)

#------------------------------------------------------------------
#assembly curve
par(mai=c(0.3,0.55,0.5,0.25))

##rank-order abundance curves for C1
C1=abundance[2:66,2:71]
timepoints=paste0(C1$Time_d,"d")
RA=C1[,3:70]
rownames(RA)=timepoints
x1=sort(RA[1,],decreasing = TRUE)

plot(x=1:68,y=x1*100,log="y",type="n",ylim=c(0.0001,120),yaxt="n",xaxt="n",xlab="",ylab="Relative cell abundance per subcommunity (%)",font.lab=2,cex.lab=1.2)
axis(1,at=c(1,68),labels=c("High","Low"))
axis(2,at=c(0.001,0.01,0.1,1,10,100),labels=c("0.001","0.01","0.1","1","10","100"),las=2)
text(x=60,y=100,labels = "C1",font = 2,cex=1.5)
abline(h=0.01,col="black",lty=1,lwd=1)
###samples picked after the adaptation phase
for(i in 7:65){
  y=sort(RA[i,],decreasing = TRUE)
  points(x=1:68,y=y*100,col="grey",type="l",lty=2)
}

###Samples picked in the adaptation phase
for(i in 1:6){
  y=sort(RA[i,],decreasing = TRUE)
  points(x=1:68,y=y*100,col=i+1,type="l",lwd=2)
}

mtext("C", side = 3, line = -41,adj=0,outer = TRUE,cex=1.5)

##rank-order abundance curves for C2
C2=abundance[67:131,2:71]
RA=C2[,3:70]
rownames(RA)=timepoints
x1=sort(RA[1,],decreasing = TRUE)
plot(x=1:68,y=x1*100,log="y",type="n",yaxt="n",ylim=c(0.0001,120),xaxt="n",xlab="",ylab="",font.lab=2)
axis(1,at=c(1,68),labels=c("High","Low"))
axis(2,at=c(0.001,0.01,0.1,1,10,100),labels=c("0.001","0.01","0.1","1","10","100"),las=2)
abline(h=0.01,col="black",lty=1,lwd=1)

text(x=60,y=100,labels = "C2",font = 2,cex=1.5)
for(i in 7:65){
  y=sort(RA[i,],decreasing = TRUE)
  points(x=1:68,y=y*100,col="grey",type="l",lty=2)
}
for(i in 1:6){
  y=sort(RA[i,],decreasing = TRUE)
  points(x=1:68,y=y*100,col=i+1,type="l",lwd=2)
}

#Legends for rank-order abundance curves
par(mai=c(0,0,0.2,0))
plot.new()
title("Rank subcommunities by relative abundance",cex.main=1.2)
legend(title="Sampling time",x="center",y="top",ncol=5,legend=c("0.25d","1d","2d","3d","5d","7d","After 7d"),col=c(2,3,4,5,6,7,"grey"),lty=c(1,1,1,1,1,1,2),lwd=1.5,text.font=2,title.adj = 0,box.col = "white",cex=1.2)
dev.off()
#Fig. 1 output as a PDF file names Figure_1.pdf

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#2.3.2 Figure 2

##rank dominant subcommunities in C1 by their abudance
RA=C1[,3:70]
rownames(RA)=timepoints
gates=colnames(RA)


for(i in 1:65){
  for(j in 1:68){
    if(RA[i,j]<0.01471){
      RA[i,j]=0
    }
  }
}

sortRA=RA
sortabundance=RA
colnames(sortRA)=1:68
colnames(sortabundance)=1:68

for(i in 1:65){
  c=sort(RA[i,],decreasing = TRUE)
  for(j in 1:68){
    if(c[j]>0){
      c[j]=names(c[j])
    }else{
      c[j]="NA"
    }
  }
  
  sortRA[i,]=c
}

for(i in 1:65){
  c=sort(RA[i,],decreasing = TRUE)
  sortabundance[i,]=c
}

##restructure the data set of C1 for plotting
data=data.frame(time=1:65,sortRA)
colnames(data)=c("time",1:68)
meltdata=melt(data,id="time")
a=as.numeric(meltdata$variable)
meltdata$variable=a


data=data.frame(time=1:65,sortabundance)
colnames(data)=c("time",1:68)
meltabundance=melt(data,id="time")
b=as.numeric(meltabundance$variable)
meltabundance$variable=b
C1melt=cbind(meltdata,abundance=meltabundance$value,reactor=rep("C1",4420))

##Rank dominant subcommunities in C2 by their abundance
RA=C2[,3:70]
rownames(RA)=paste0(timepoints,"d")
gates=colnames(RA)

#
for(i in 1:65){
  for(j in 1:68){
    if(RA[i,j]<0.01471){
      RA[i,j]=0
    }
  }
}

sortRA=RA
sortabundance=RA
colnames(sortRA)=1:68
colnames(sortabundance)=1:68

for(i in 1:65){
  c=sort(RA[i,],decreasing = TRUE)
  for(j in 1:68){
    if(c[j]>0){
      c[j]=names(c[j])
    }else{
      c[j]="NA"
    }
  }
  
  sortRA[i,]=c
}

for(i in 1:65){
  c=sort(RA[i,],decreasing = TRUE)
  sortabundance[i,]=c
}

data=data.frame(time=1:65,sortRA)
colnames(data)=c("time",1:68)
meltdata=melt(data,id="time")
a=as.numeric(meltdata$variable)
meltdata$variable=a
##restructure the data set of C2 for plotting
data=data.frame(time=1:65,sortabundance)
colnames(data)=c("time",1:68)
meltabundance=melt(data,id="time")
b=as.numeric(meltabundance$variable)
meltabundance$variable=b
C2melt=cbind(meltdata,abundance=meltabundance$value,reactor=rep("C2",4420))

##cobined data sets of C1melt and C2melt and only keep dominant subcommunities for visulization
bindC1C2=rbind(C1melt,C2melt)
meltdata=bindC1C2[1,]
for(i in 2:8840){
  if(bindC1C2$value[i]!="NA"){
    meltdata=rbind(meltdata,bindC1C2[i,])
  }
}
c=meltdata$abundance*100#convert values to persentage e.g. 0.01=1%
meltdata$abundance=c
reactorname=c('C1'="C1",'C2'="C2")
bre=c(9,36,64)
textsize=15
axisize=15
meltdata=meltdata[order(meltdata$time),]

##visulzation with ggplot2
### Fig. 2A
colorcode=rainbow(68)
p=ggplot(data=meltdata,aes(x=time,y=abundance,fill=value,group=variable))+geom_bar(position = position_stack(reverse = TRUE),colour="black",alpha=0.7,stat="identity")+facet_grid(.~reactor,labeller = as_labeller(reactorname))+scale_fill_manual(values=colorcode,name="Color key for subcommunities")
q=p+scale_x_continuous(limits = c(0,67),expand = c(0,0.025),breaks=c(bre),labels=c(timepoints[bre]))+scale_y_continuous(limits=c(0,100),expand=c(0,0.025),name = "Relative abundance (%)")+theme(legend.position = "top",legend.direction = "vertical") + guides(fill = guide_legend(ncol = 34))
r=q+theme(legend.text = element_text(size=12),legend.title = element_text(face="bold",size=textsize),legend.key.size = unit(0.35,"cm"),panel.background = element_blank())+xlab("Time (d)")+theme(text=element_text(size=axisize),axis.text.x = element_text(color = "black",hjust=0.5,vjust = 0.5))+theme(axis.text.y = element_text(color = "black"))+theme(strip.background =element_rect(colour = "grey"),panel.border = element_rect(color = "grey",fill = NA),strip.text = element_text(face="bold",size=textsize),axis.title = element_text(face="bold",size=textsize))
r
###Fig. 2B
div=divCs
a=ggplot(data=div,aes(color=reactor))+facet_grid(.~reactor,labeller = as_labeller(reactorname))+geom_line(data=div[259:260,],aes(x=time,y=Values),linetype=4)+geom_line(data=div[261:262,],aes(x=time,y=Values),linetype=4)+geom_line(data=div[1:65,],aes(x=time,y=Values),size=1)+geom_line(data=div[66:129,],aes(x=time,y=Values),size=1,linetype=2)+geom_point(data=div[1:129,],aes(x=time,y=Values,shape=group),fill="white",size=2.5,stroke =1.2)+geom_point(data=div[130:194,],aes(x=time,y=Values),size=3)+geom_line(data=div[130:194,],aes(x=time,y=Values),size=1)+geom_line(data=div[195:258,],aes(x=time,y=Values),size=1,linetype=2)+geom_point(data=div[195:258,],aes(x=time,y=Values),fill="white",shape=21,size=2.5,stroke=1.2)
b=a+scale_color_manual(values = c("gray50","grey"))+scale_x_continuous(limits = c(0,67),expand = c(0,0.025),breaks=c(bre),labels=c(timepoints[bre]))+scale_y_continuous(limits=c(0,30),expand=c(0.03,0.025),breaks=c(0,5,10,15,20,25,30),labels=c(0,5,10,15,20,25,30),name = "Diversity values")+xlab("Time (d)")+theme(text = element_text(size=axisize),axis.text.x = element_text(color="black",hjust=0.5,vjust = 0.5))+theme(axis.text.y = element_text(color = "black"))+theme(strip.text.x =element_blank(),strip.background =element_blank(),panel.border = element_rect(color = "grey",fill = NA),axis.title = element_text(face="bold",size=textsize))+theme(panel.background =element_blank(),plot.margin =unit(c(0.5,0.15,0.2,0.5),"cm"),plot.background = element_rect(fill = "NA"))+scale_shape_manual(values=c(19,21))+theme(legend.position = "none")
b

#output as a PDF file which named as Figure_2.pdf
pdf("Figure_2.pdf",width =15,height=15)
grid.newpage()
pushViewport(viewport(layout=grid.layout(51,1)))
print(r, vp=viewport(layout.pos.row=c(2:26),layout.pos.col=1))
print(b, vp=viewport(layout.pos.row=c(27:51),layout.pos.col=1))
grid.text(expression(bold("A")),x=unit(0.02,"npc"),y=unit(0.99,"npc"),gp=gpar(fontsize=20))
grid.text(expression(bold("B")),x=unit(0.02,"npc"),y=unit(0.48,"npc"),gp=gpar(fontsize=20))
dev.off()
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
#2.3.3 Figure 3

#Command for creating pdf file
pdf("Figure_3.pdf",height=9.8,width=10.1)
layout(matrix(c(1,1,1,2,2,2,3,3,4,4,5,5,6,6,6,6,6,6,6,6,6,6,6,6), ncol=6, byrow=TRUE), heights=c(7.5,5.1,0.9))
par(mai=c(0.55,0.55,0.25,0.25),omi=c(0,0.1,0.3,0.1),xpd=FALSE)
#------------------------------------------------------------------
#Virsulization of the reference space
plot(x=EcoSat$time,y=EcoSat$D1,type="n",pch=19,col="gray50",ylim=c(0,1),xaxt="n",xlab="Time (d)",ylab="Deviation from reference state",font.lab=2,cex.lab=1.2)

abline(h=0.251,col="deepskyblue2",lty=2,lwd=2)
abline(h=0.59,col="brown3",lty=2,lwd=2)
abline(h=0.471,col="darkseagreen3",lty=2,lwd=2)
points(x=EcoSat$time,y=EcoSat$D1,type="o",pch=19,col="deepskyblue2",cex=1.5,lwd=2)
points(x=EcoSat$time,y=EcoSat$D2,type="o",pch=19,col="brown3",cex=1.5,lwd=2)
points(x=EcoSat$time,y=EcoSat$D3,type="o",pch=19,col="darkseagreen3",cex=1.5,lwd=2)
axis(side=1,at= seq(0,91,7))
mtext("A", side = 3, line = -1,adj=0,outer = TRUE,cex=1.5)
#------------------------------------------------------------------
#Dissimilarity analysis with relative cells abundance per subcommunities (all five reactor)
#Dissimilarity analysis performed with data from controls (see 2.1.1)

# draw nMDS
plot(mds.out, type="n",font.lab=2,cex.lab=1.2,bg="transparent")
mtext("B", side = 3, line = -1,adj=0.5,outer = TRUE,cex=1.5)
# mark reference state
arrows(mds.out$points[1,1],mds.out$points[1,2],mds.out$points[132,1],mds.out$points[132,2], length = 0.08, angle = 15, lwd=1.5, col ="deepskyblue2")
arrows(mds.out$points[1,1],mds.out$points[1,2],mds.out$points[197,1],mds.out$points[197,2], length = 0.08, angle = 15, lwd=1.5, col ="brown3")
arrows(mds.out$points[1,1],mds.out$points[1,2],mds.out$points[262,1],mds.out$points[262,2], length = 0.08, angle = 15, lwd=1.5, col ="darkseagreen3")
points(mds.out$points[1,1],mds.out$points[1,2], col ="black", pch = 19, cex = 1.5,lwd=3)


points(mds.out$points[197:261,], col ="brown3",type="l",lty=2)
points(mds.out$points[197:261,], col ="brown3", pch = 21,bg="white",cex = size,type="p",lwd=2)

points(mds.out$points[262:326,], col ="darkseagreen3",type="l",lty=2)
points(mds.out$points[262:326,], col ="darkseagreen3", pch = 21,bg="white",cex = size,type="p",lwd=2)

points(mds.out$points[132:196,], col ="deepskyblue2",  type="l",lty=2)
points(mds.out$points[132:196,], col ="deepskyblue2", pch = 21,bg="white", cex = size,type="p",lwd=2)
#------------------------------------------------------------------
#assembly curve
par(mai=c(0.4,0.55,0.3,0.05))

##Rank-order abundance curves of D1
D1=abundance[132:196,2:71]
RA=D1[,3:70]
rownames(RA)=timepoints
x1=sort(RA[1,],decreasing = TRUE)

plot(x=1:68,y=x1*100,log="y",type="n",ylim=c(0.0001,120),xaxt="n",yaxt="n",xlab="",ylab="Relative abundance per subcommunity(%)  ",font.lab=2,cex.lab=1.2)
axis(1,at=c(1,68),labels=c("High","Low"))
axis(2,at=c(0.001,0.01,0.1,1,10,100),labels=c("0.001","0.01","0.1","1","10","100"),las=2)
text(x=60,y=100,labels = "D1",font = 2,cex=1.5)
abline(h=0.01,col="black",lty=1,lwd=1)
for(i in 7:65){
  y=sort(RA[i,],decreasing = TRUE)
  points(x=1:68,y=y*100,col="grey",type="l",lty=2)
}


for(i in 1:6){
  y=sort(RA[i,],decreasing = TRUE)
  points(x=1:68,y=y*100,col=i+1,type="l",lwd=2)
}
mtext("C", side = 3, line = -39,adj=0,outer = TRUE,cex=1.5)

##Rank-order abundance curves of D2
par(mai=c(0.4,0.45,0.3,0.15))
D2=abundance[197:261,2:71]
RA=D2[,3:70]
rownames(RA)=timepoints
x1=sort(RA[1,],decreasing = TRUE)

plot(x=1:68,y=x1*100,log="y",type="n",ylim=c(0.0001,120),xaxt="n",yaxt="n",xlab="",ylab=" ",font.lab=2)
axis(1,at=c(1,68),labels=c("High","Low"))
axis(2,at=c(0.001,0.01,0.1,1,10,100),labels=c("0.001","0.01","0.1","1","10","100"),las=2)
text(x=60,y=100,labels = "D2",font = 2,cex=1.5)
abline(h=0.01,col="black",lty=1,lwd=1)
for(i in 7:65){
  y=sort(RA[i,],decreasing = TRUE)
  points(x=1:68,y=y*100,col="grey",type="l",lty=2)
}


for(i in 1:6){
  y=sort(RA[i,],decreasing = TRUE)
  points(x=1:68,y=y*100,col=i+1,type="l",lwd=2)
}

##Rank-order abundance curves of D3
par(mai=c(0.4,0.35,0.3,0.25))
D3=abundance[262:326,2:71]
RA=D3[,3:70]
rownames(RA)=timepoints
x1=sort(RA[1,],decreasing = TRUE)

plot(x=1:68,y=x1*100,log="y",type="n",ylim=c(0.0001,120),xaxt="n",yaxt="n",xlab="",ylab=" ",font.lab=2)
axis(1,at=c(1,68),labels=c("High","Low"))
axis(2,at=c(0.001,0.01,0.1,1,10,100),labels=c("0.001","0.01","0.1","1","10","100"),las=2)
text(x=60,y=100,labels = "D3",font = 2,cex=1.5)
abline(h=0.01,col="black",lty=1,lwd=1)
for(i in 7:65){
  y=sort(RA[i,],decreasing = TRUE)
  points(x=1:68,y=y*100,col="grey",type="l",lty=2)
}


for(i in 1:6){
  y=sort(RA[i,],decreasing = TRUE)
  points(x=1:68,y=y*100,col=i+1,type="l",lwd=2)
}

#Legends for rank-order abundance curvespar(mai=c(0.3,0,0.2,0))
plot.new()
title("Rank subcommunities by relative abundance",cex.main=1.2)
legend(title="Sampling time",x="center",y="top",ncol=5,legend=c("0.25d","1d","2d","3d","5d","7d","After 7d"),col=c(2,3,4,5,6,7,"grey"),lty=c(1,1,1,1,1,1,2),lwd=1.5,text.font=2,title.adj = 0,box.col = "white",cex=1.2)
dev.off()
#Fig. 3 output as a PDF file names Figure_3.pdf
#-------------------------------------------------------------
#-------------------------------------------------------------
#2.3.4 Figure 4
##Rank dominant subcommunities in D1 by their abundances
RA=D1[,3:70]
rownames(RA)=timepoints
gates=colnames(RA)

for(i in 1:65){
  for(j in 1:68){
    if(RA[i,j]<0.01471){
      RA[i,j]=0
    }
  }
}

sortRA=RA
sortabundance=RA
colnames(sortRA)=1:68
colnames(sortabundance)=1:68

for(i in 1:65){
  c=sort(RA[i,],decreasing = TRUE)
  for(j in 1:68){
    if(c[j]>0){
      c[j]=names(c[j])
    }else{
      c[j]="NA"
    }
  }
  
  sortRA[i,]=c
}

for(i in 1:65){
  c=sort(RA[i,],decreasing = TRUE)
  sortabundance[i,]=c
}


data=data.frame(time=1:65,sortRA)
colnames(data)=c("time",1:68)
meltdata=melt(data,id="time")
a=as.numeric(meltdata$variable)
meltdata$variable=a

##restructure the data set of D1 for plotting
data=data.frame(time=1:65,sortabundance)
colnames(data)=c("time",1:68)
meltabundance=melt(data,id="time")
b=as.numeric(meltabundance$variable)
meltabundance$variable=b
D1melt=cbind(meltdata,abundance=meltabundance$value,reactor=rep("D1",4420))

#Rank dominant subcommunities in D2 by their abundances
RA=D2[,3:70]
rownames(RA)=timepoints
gates=colnames(RA)

for(i in 1:65){
  for(j in 1:68){
    if(RA[i,j]<0.01471){
      RA[i,j]=0
    }
  }
}

sortRA=RA
sortabundance=RA
colnames(sortRA)=1:68
colnames(sortabundance)=1:68

for(i in 1:65){
  c=sort(RA[i,],decreasing = TRUE)
  for(j in 1:68){
    if(c[j]>0){
      c[j]=names(c[j])
    }else{
      c[j]="NA"
    }
  }
  
  sortRA[i,]=c
}

for(i in 1:65){
  c=sort(RA[i,],decreasing = TRUE)
  sortabundance[i,]=c
}

data=data.frame(time=1:65,sortRA)
colnames(data)=c("time",1:68)
meltdata=melt(data,id="time")
a=as.numeric(meltdata$variable)
meltdata$variable=a

##restructure the data set of D2 for plotting
data=data.frame(time=1:65,sortabundance)
colnames(data)=c("time",1:68)
meltabundance=melt(data,id="time")
b=as.numeric(meltabundance$variable)
meltabundance$variable=b
D2melt=cbind(meltdata,abundance=meltabundance$value,reactor=rep("D2",4420))

##Rank dominant subcommunities in D3 by their abundances
RA=D3[,3:70]
rownames(RA)=timepoints
gates=colnames(RA)

for(i in 1:65){
  for(j in 1:68){
    if(RA[i,j]<0.01471){
      RA[i,j]=0
    }
  }
}

sortRA=RA
sortabundance=RA
colnames(sortRA)=1:68
colnames(sortabundance)=1:68

for(i in 1:65){
  c=sort(RA[i,],decreasing = TRUE)
  for(j in 1:68){
    if(c[j]>0){
      c[j]=names(c[j])
    }else{
      c[j]="NA"
    }
  }
  
  sortRA[i,]=c
}

for(i in 1:65){
  c=sort(RA[i,],decreasing = TRUE)
  sortabundance[i,]=c
}


data=data.frame(time=1:65,sortRA)
colnames(data)=c("time",1:68)
meltdata=melt(data,id="time")
a=as.numeric(meltdata$variable)
meltdata$variable=a

##restructure the data set of D3 for plotting
data=data.frame(time=1:65,sortabundance)
colnames(data)=c("time",1:68)
meltabundance=melt(data,id="time")
b=as.numeric(meltabundance$variable)
meltabundance$variable=b
D3melt=cbind(meltdata,abundance=meltabundance$value,reactor=rep("D3",4420))

##cobined data sets of C1melt and C2melt and only keep dominant subcommunities for visulization
bindD1D2D3=rbind(D1melt,D2melt,D3melt)
meltdata=bindD1D2D3[1,]
for(i in 2:13260){
  if(bindD1D2D3$value[i]!="NA"){
    meltdata=rbind(meltdata,bindD1D2D3[i,])
  }
}
c=meltdata$abundance*100
meltdata$abundance=c

##visulzation with ggplot2
### Fig. 4A
reactorname=c('D1'="D1",'D2'="D2",'D3'="D3")
graphname=c('1'="Alpha and inter-beta diversity",'2'="Intra-beta diversity")

bre=c(9,36,64)
textsize=18

p=ggplot(data=meltdata,aes(x=time,y=abundance,fill=value,group=variable))+geom_bar(position = position_stack(reverse = TRUE),colour="black",alpha=0.70,stat="identity")+facet_grid(.~reactor,labeller = as_labeller(reactorname))+scale_fill_manual(values=colorcode,name="Color key for subcommunities")
q=p+scale_x_continuous(limits = c(0,67),expand = c(0,0.025),breaks=c(bre),labels=c(timepoints[bre]))+scale_y_continuous(limits=c(0,100),expand=c(0,0.025),name = "Relative abundance (%)")+theme(legend.position = "top",legend.direction = "vertical") + guides(fill = guide_legend(ncol = 34))
r=q+theme(legend.text = element_text(size=12.5),legend.title = element_text(face="bold",size=textsize),legend.key.size = unit(0.35,"cm"),panel.background = element_blank())+xlab("Time (d)")+theme(text=element_text(size=textsize),axis.text.x = element_text(color = "black",hjust=0.5,vjust = 0.5))+theme(axis.text.y = element_text(color = "black"))+theme(strip.background =element_rect(colour = "grey"),panel.border = element_rect(color = "grey",fill = NA),strip.text = element_text(face="bold",size=textsize),axis.title = element_text(face="bold",size=textsize))
r

### Fig. 4B
div=divDs
a=ggplot(data=div[1:387,],aes(shape=as.factor(group),color=reactor,linetype=as.factor(group)))+facet_grid(.~reactor,labeller = as_labeller(reactorname))+geom_hline(yintercept = 8,color="gray50",size=1,linetype=4)+geom_line(data=div[1:387,],aes(x=time,y=Values),size=1)+scale_linetype_manual(values=c(1,2))+scale_color_manual(values = c("#00b2ee","#cd3333","#9bcd9b"))+scale_shape_manual(values=c(19,21))+geom_point(aes(x=time,y=Values),fill="white",size=2.5,stroke =1.2)+scale_x_continuous(limits = c(0,67),expand = c(0,0.025),breaks=c(bre),labels=c(timepoints[bre]))+scale_y_continuous(limits=c(0,30),expand=c(0.03,0.025),breaks=c(0,5,10,15,20,25,30),labels=c(0,5,10,15,20,25,30),name = "Diversity values")+xlab("Time (d)")
b=a+theme(plot.margin =unit(c(0.3,0.15,0.2,0.7),"cm"),legend.position = "none",strip.background = element_blank(),strip.text = element_blank())+theme(text = element_text(size=textsize),axis.text.x = element_text(hjust=0.5,vjust = 0.5,color = "black"),axis.text.y = element_text(color = "black"),panel.border = element_rect(color = "grey",fill = NA),axis.title = element_text(face="bold"),panel.background =element_blank(),plot.background = element_rect(fill = "NA"))
b

###Temperature marker
labels=c(expression(30~degree*C),expression(40~degree*C))
s=ggplot(data=te)+facet_grid(.~reactor)+geom_line(data=te,aes(x=time,y=Values),linetype=1,size=1)+scale_x_continuous(limits = c(0,67),expand = c(0,0.025),breaks=c(bre),labels=c(timepoints[bre]))+scale_y_continuous(expand=c(0,0),breaks=c(30,40),labels=labels)+theme_classic()
t=s+theme(plot.margin =unit(c(0.2,0.15,0.15,0.66),"cm"),axis.line = element_blank(),axis.ticks = element_blank(),strip.background = element_blank(),strip.text = element_blank(),axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=15,color = "black"))
t
#Output Fig. 4 as a PDF file named Figure_4.pdf
pdf("Figure_4.pdf",width =22.5,height=15.75)
grid.newpage()
pushViewport(viewport(layout=grid.layout(21,1)))
print(r, vp=viewport(layout.pos.row=c(1:10),layout.pos.col=1))
print(t, vp=viewport(layout.pos.row=11,layout.pos.col=1))
print(b, vp=viewport(layout.pos.row=c(12:21),layout.pos.col=1))
grid.text(expression(bold("A")),x=unit(0.01,"npc"),y=unit(0.99,"npc"),gp=gpar(fontsize=22))
grid.text(expression(bold("B")),x=unit(0.01,"npc"),y=unit(0.46,"npc"),gp=gpar(fontsize=22))
dev.off()
#-----------------------------------------------------------------


#============================================================================
#2.4 Codes for the supporting information
#2.4.1 S3 Abiotic parameters
mydata=abio
melt=melt(mydata,id=c("Group","Time"))

colnames(melt)=c("group","x","parameter","y")
melt=na.omit(melt)
Nrow=nrow(melt)

melt_tem=melt[1:390,]
melt_abio=melt[391:Nrow,]

labelabio=c(expression(paste("T (",degree,"C)")),expression("pH"),expression(paste("EC (mS cm"^"-1",")")),expression(paste("NH4 (mg N L"^"-1",")")),expression(paste("PHO (mg P L"^"-1",")")),expression(paste("CODs (mg L"^"-1",")")),expression(paste("CODt (mg L"^"-1",")")),expression(paste("CODb (mg L"^"-1",")")),expression(paste("Number of cells ", "(10"^"6"," mL"^"-1",")")),expression("OD"))
p=ggplot(data=melt,aes(x=x,y=y,colour=parameter))+facet_grid(parameter~group,scales = "free_y")+geom_point(data=melt_tem,size=1.25)+geom_point(data=melt_abio,size=1.25)+geom_line(data=melt_tem,size=0.75)+geom_line(data=melt_abio,size=0.75)+scale_color_discrete(name="Parameters",breaks=colnames(mydata)[3:ncol(mydata)],labels=labelabio)
q=p+labs(x="Time (d)",y="Measurement of biotic and abiotic parameters")+theme(strip.text = element_text(size=15,face = "bold"),text = element_text(face="plain"),axis.title = element_text(size=15,face = "bold"),axis.text = element_text(size=10),legend.title = element_text(size=13,face = "bold"),legend.text = element_text(size=13,face = "plain"),legend.text.align = 0,strip.text.y = element_blank(),legend.position = "bottom")
q
#Output as a PDF file named AbioParameters.pdf
pdf("AbioParameters.pdf",height = 12,width = 12)
print(q)
dev.off()
#-----------------------------------------------------------------
#2.4.2 Analysis of dis-/similarity (ANOSIM)
##ANOSIM for the whole data set of C1, C2, D1, D2 and D3.
###subsampled in groups
bcdist=vegdist(abundance[2:326,4:71],"bray")
anosim_bc_group=anosim(bcdist,group[2:326])
plot(anosim_bc_group)
summary(anosim_bc_group)
###subsampled in reactors
anosim_bc_reactor=anosim(bcdist,reactor[2:326])
plot(anosim_bc_reactor)
summary(anosim_bc_reactor)

##ANOSIM for pairwise reactors in the adaptation phase.(Tab. S8.1)
#C1C2
dist_inADP_C1C2=vegdist(abundance[c(2:7,67:72),4:71],"bray")
anosim_inADP_C1C2=anosim(dist_inADP_C1C2,reactor[c(2:7,67:72)])
plot(anosim_inADP_C1C2)
summary(anosim_inADP_C1C2)
#C1D1
dist_inADP_C1D1=vegdist(abundance[c(2:7,132:137),4:71],"bray")
anosim_inADP_C1D1=anosim(dist_inADP_C1D1,reactor[c(2:7,132:137)])
plot(anosim_inADP_C1D1)
summary(anosim_inADP_C1D1)
#C1D2
dist_inADP_C1D2=vegdist(abundance[c(2:7,197:202),4:71],"bray")
anosim_inADP_C1D2=anosim(dist_inADP_C1D2,reactor[c(2:7,197:202)])
plot(anosim_inADP_C1D2)
summary(anosim_inADP_C1D2)
#C1D3
dist_inADP_C1D3=vegdist(abundance[c(2:7,262:267),4:71],"bray")
anosim_inADP_C1D3=anosim(dist_inADP_C1D3,reactor[c(2:7,262:267)])
plot(anosim_inADP_C1D3)
summary(anosim_inADP_C1D3)
#C2D1
dist_inADP_C2D1=vegdist(abundance[c(67:72,132:137),4:71],"bray")
anosim_inADP_C2D1=anosim(dist_inADP_C2D1,reactor[c(67:72,132:137)])
plot(anosim_inADP_C2D1)
summary(anosim_inADP_C2D1)
#C2D2
dist_inADP_C2D2=vegdist(abundance[c(67:72,197:202),4:71],"bray")
anosim_inADP_C2D2=anosim(dist_inADP_C2D2,reactor[c(67:72,197:202)])
plot(anosim_inADP_C2D2)
summary(anosim_inADP_C2D2)
#C2D3
dist_inADP_C2D3=vegdist(abundance[c(67:72,262:267),4:71],"bray")
anosim_inADP_C2D3=anosim(dist_inADP_C2D3,reactor[c(67:72,262:267)])
plot(anosim_inADP_C2D3)
summary(anosim_inADP_C2D3)
#D1D2
dist_inADP_D1D2=vegdist(abundance[c(132:137,197:202),4:71],"bray")
anosim_inADP_D1D2=anosim(dist_inADP_D1D2,reactor[c(132:137,197:202)])
plot(anosim_inADP_D1D2)
summary(anosim_inADP_D1D2)
#D1D3
dist_inADP_D1D3=vegdist(abundance[c(132:137,262:267),4:71],"bray")
anosim_inADP_D1D3=anosim(dist_inADP_D1D3,reactor[c(132:137,262:267)])
plot(anosim_inADP_D1D3)
summary(anosim_inADP_D1D3)
#D2D3
dist_inADP_D2D3=vegdist(abundance[c(197:202,262:267),4:71],"bray")
anosim_inADP_D2D3=anosim(dist_inADP_D2D3,reactor[c(197:202,262:267)])
plot(anosim_inADP_D2D3)
summary(anosim_inADP_D2D3)

##ANOSIM for pairwise reactors after the adaptation phase.(Tab. S8.1)
#C1C2
dist_afterADP_C1C2=vegdist(abundance[c(8:66,73:131),4:71],"bray")
anosim_afterADP_C1C2=anosim(dist_afterADP_C1C2,reactor[c(8:66,73:131)])
plot(anosim_afterADP_C1C2)
summary(anosim_afterADP_C1C2)
#C1D1
dist_afterADP_C1D1=vegdist(abundance[c(8:66,138:196),4:71],"bray")
anosim_afterADP_C1D1=anosim(dist_afterADP_C1D1,reactor[c(8:66,138:196)])
plot(anosim_afterADP_C1D1)
summary(anosim_afterADP_C1D1)
#C1D2
dist_afterADP_C1D2=vegdist(abundance[c(8:66,203:261),4:71],"bray")
anosim_afterADP_C1D2=anosim(dist_afterADP_C1D2,reactor[c(8:66,203:261)])
plot(anosim_afterADP_C1D2)
summary(anosim_afterADP_C1D2)
#C1D3
dist_afterADP_C1D3=vegdist(abundance[c(8:66,268:326),4:71],"bray")
anosim_afterADP_C1D3=anosim(dist_afterADP_C1D3,reactor[c(8:66,268:326)])
plot(anosim_afterADP_C1D3)
summary(anosim_afterADP_C1D3)
#C2D1
dist_afterADP_C2D1=vegdist(abundance[c(73:131,138:196),4:71],"bray")
anosim_afterADP_C2D1=anosim(dist_afterADP_C2D1,reactor[c(73:131,138:196)])
plot(anosim_afterADP_C2D1)
summary(anosim_afterADP_C2D1)
#C2D2
dist_afterADP_C2D2=vegdist(abundance[c(73:131,203:261),4:71],"bray")
anosim_afterADP_C2D2=anosim(dist_afterADP_C2D2,reactor[c(73:131,203:261)])
plot(anosim_afterADP_C2D2)
summary(anosim_afterADP_C2D2)
#C2D3
dist_afterADP_C2D3=vegdist(abundance[c(73:131,268:326),4:71],"bray")
anosim_afterADP_C2D3=anosim(dist_afterADP_C2D3,reactor[c(73:131,268:326)])
plot(anosim_afterADP_C2D3)
summary(anosim_afterADP_C2D3)
#D1D2
dist_afterADP_D1D2=vegdist(abundance[c(138:196,203:261),4:71],"bray")
anosim_afterADP_D1D2=anosim(dist_afterADP_D1D2,reactor[c(138:196,203:261)])
plot(anosim_afterADP_D1D2)
summary(anosim_afterADP_D1D2)
#D1D3
dist_afterADP_D1D3=vegdist(abundance[c(138:196,268:326),4:71],"bray")
anosim_afterADP_D1D3=anosim(dist_afterADP_D1D3,reactor[c(138:196,268:326)])
plot(anosim_afterADP_D1D3)
summary(anosim_afterADP_D1D3)
#D2D3
dist_afterADP_D2D3=vegdist(abundance[c(203:261,268:326),4:71],"bray")
anosim_afterADP_D2D3=anosim(dist_afterADP_D2D3,reactor[c(203:261,268:326)])
plot(anosim_afterADP_D2D3)
summary(anosim_afterADP_D2D3)
#-----------------------------------------------------------------
#2.4.3 Correlation analysis
library(Hmisc)        
library(psych)   
library(gplots)              
library(qgraph)


Bll=All
#only keep dominant subcommunities
for(i in 1:325){
  for(j in 2:69){
    if(All[i,j]<=0.0147){
      Bll[i,j]=0
    }
  }
}

#function used for visulization correlation network
#s:the No.range of subsampled samples, e.g. for the first phase per reactor, s=1:10.
#parameters: the numbers of columns which are parameters used for correlation anslysis, e.g. all parameters, parameters=2:79
#min: the threshold for determing significant correlations, while in this study min=0.75, that means only correlations values above 0.75 or below -0.75 were determined as dominant correlations.
#This function output correlations as a network include given samples and also report the numbers of significant correlations, while 'BB':Sc vs. Sc; 'AB': SC vs. Abio; and 'AA': Abio vs. Abio.
#calculation based on dominant subcommunities only.
min=0.75
CP=function(s,parameters,min){
  col=c(rep("deepskyblue3",68),rep("brown3",length(parameters)-68))
  allmatrix<-do.call(cbind,Bll[s,parameters])
  j=0
  for(i in 1:68){
    if(sum(allmatrix[,i])<=0.0147){
      j=c(j,i)
    }
  }
  j=c(j[-1])
  bionumber=68-length(j)
  allmatrix=allmatrix[,-j]
  col=col[-j]
  n<-colnames(allmatrix)
  
  
  
  cor<-rcorr(x=allmatrix,type=c("spearman"))
  r<-cor$r
  p<-c(cor$P)
  colnames(r)<-n
  rownames(r)<-n
  bh<-sapply(as.numeric(p), p.adjust, method="BH")
  padj<-c(bh)
  padj[is.na(padj)]<-0
  
  cor1=rcorr(x=allmatrix,type=c("spearman"))
  r1<-cor1$r
  p1<-c(cor1$P)
  colnames(r1)<-n
  rownames(r1)<-n
  bh1<-sapply(as.numeric(p1), p.adjust, method="BH")
  padj1<-c(bh1)
  padj1[is.na(padj1)]<-0
  
  rr1=r
  for(i in 1:length(n)){
    for(j in 1:length(n)){
      if(r[i,j]=="NaN"|r1[i,j]=="NaN"){
        rr1[i,j]=0
      }else{
        if(i>j){
          if(r[i,j]>=min&&r1[i,j]>=min){
            rr1[i,j]=min
          }else{
            if(r[i,j]<=-min&&r1[i,j]<=-min){
              rr1[i,j]=-min
            }else{
              rr1[i,j]=0
            }
          }
        }else{rr1[i,j]=0}
      }
    }
  }
  for (i in 1:length(padj)) {
    if (padj[i]>0.05){rr1[[i]]<-0}
  }
  
  for (i in 1:length(padj1)) {
    if (padj1[i]>0.05){rr1[[i]]<-0}
  }
  BB=0
  for(i in 1:bionumber){
    for(j in 1:bionumber){
      if(i!=j){
        if(rr1[i,j]>=min|rr1[i,j]<=-min){
          BB=BB+1
        }
      }
    }
  }
  AB=0
  for(i in (bionumber+1):length(n)){
    for(j in 1:bionumber){
      if(rr1[i,j]>=min|rr1[i,j]<=-min)
      {
        AB=AB+1
      }
    }
  }
  AA=0
  for(i in (bionumber+1):length(n)){
    for(j in (bionumber+1):length(n)){
      if(rr1[i,j]>=min|rr1[i,j]<=-min)
      {
        AA=AA+1
      }
    }
  }
  
  p=qgraph(rr1,vsize=3.5,legend=FALSE,borders=TRUE,layout="spring",parallelEdge=TRUE,details=FALSE,negCol="#91cf60",posCol="#ef8a62",esize=1,directed=FALSE,color=col,label.scale=TRUE,label.font=2,mar=c(1,3,3,1))
  out=c('BB'=BB,'AB'=AB,'AA'=AA)
  return(out)
}
par(mfrow=c(1,1),omi=c(0,0,0,0))

#Significant correlations were visulized as net works
## 1st phase in C1 
tiff("C1_1st.tif",height=3,width=3,units = "in",res = 600,compression = "lzw")
s=c(1:10)
parameters=c(2:78)
counts=CP(s,parameters,min)
cornum=data.frame('C1P1'=counts)
dev.off()
## 2nd phase drift 1 in C1
tiff("C1_drift1.tif",height=3,width=3,units = "in",res = 600,compression = "lzw")
s=c(12:21)
parameters=c(2:78)
counts=CP(s,parameters,min)
cornum$C1P2=counts
dev.off()
## 2nd phase drift 2 in C1
tiff("C1_drift2.tif",height=3,width=3,units = "in",res = 600,compression = "lzw")
s=c(34:43)
parameters=c(2:78)
counts=CP(s,parameters,min)
cornum$C1P3=counts
dev.off()
## 3rd phase in C1
tiff("C1_last.tif",height=3,width=3,units = "in",res = 600,compression = "lzw")
s=c(56:65)
parameters=c(2:78)
counts=CP(s,parameters,min)
cornum$C1P6=counts
dev.off()
#----------------------------------------
## 1st phase in C2
tiff("C2_1st.tif",height=3,width=3,units = "in",res = 600,compression = "lzw")
s=c(1:10)+65
parameters=c(2:78)
counts=CP(s,parameters,min)
cornum$C2P1=counts
dev.off()
## 2nd phase drift 1 in C2
tiff("C2_drift1.tif",height=3,width=3,units = "in",res = 600,compression = "lzw")
s=c(8:17)+65
parameters=c(2:78)
counts=CP(s,parameters,min)
cornum$C2P2=counts
dev.off()
## 2nd phase drift 2 in C2
tiff("C2_drift2.tif",height=3,width=3,units = "in",res = 600,compression = "lzw")
s=c(35:44)+65
parameters=c(2:78)
counts=CP(s,parameters,min)
cornum$C2P3=counts
dev.off()
## 2nd phase drift 3 in C2
tiff("C2_drift3.tif",height=3,width=3,units = "in",res = 600,compression = "lzw")
s=c(45:54)+65
parameters=c(2:78)
counts=CP(s,parameters,min)
cornum$C2P4=counts
dev.off()
## 3rd phase in C2
tiff("C2_last.tif",height=3,width=3,units = "in",res = 600,compression = "lzw")
s=c(56:65)+65
parameters=c(2:78)
counts=CP(s,parameters,min)
cornum$C2P6=counts
dev.off()
#----------------------------------------
## 1st phase in D1
tiff("D1_1st.tif",height=3,width=3,units = "in",res = 600,compression = "lzw")
s=c(1:10)+130
parameters=c(2:78)
counts=CP(s,parameters,min)
cornum$D1P1=counts
dev.off()
## 2nd phase T1 in D1
tiff("D1_T1.tif",height=3,width=3,units = "in",res = 600,compression = "lzw")
s=c(6:15)+130
parameters=c(2:78)
counts=CP(s,parameters,min)
cornum$D1P2=counts
dev.off()
## 2nd phase T2 in D1
tiff("D1_T2.tif",height=3,width=3,units = "in",res = 600,compression = "lzw")
s=c(22:31)+130
parameters=c(2:78)
counts=CP(s,parameters,min)
cornum$D1P3=counts
dev.off()
## 2nd phase T3 in D1
tiff("D1_T3.tif",height=3,width=3,units = "in",res = 600,compression = "lzw")
s=c(35:44)+130
parameters=c(2:78)
counts=CP(s,parameters,min)
cornum$D1P4=counts
dev.off()
## 2nd phase T4 in D1
tiff("D1_T4.tif",height=3,width=3,units = "in",res = 600,compression = "lzw")
s=c(49:58)+130
parameters=c(2:78)
counts=CP(s,parameters,min)
cornum$D1P5=counts
dev.off()
## 3rd phase in D1
tiff("D1_last.tif",height=3,width=3,units = "in",res = 600,compression = "lzw")
s=c(56:65)+130
parameters=c(2:78)
counts=CP(s,parameters,min)
cornum$D1P6=counts
dev.off()
#----------------------------------------
## 1st phase in D2
tiff("D2_1st.tif",height=3,width=3,units = "in",res = 600,compression = "lzw")
s=c(1:10)+195
parameters=c(2:78)
counts=CP(s,parameters,min)
cornum$D2P1=counts
dev.off()
## 2nd phase T1 in D2
tiff("D2_T1.tif",height=3,width=3,units = "in",res = 600,compression = "lzw")
s=c(6:15)+195
parameters=c(2:78)
counts=CP(s,parameters,min)
cornum$D2P2=counts
dev.off()
## 2nd phase T2 in D2
tiff("D2_T2.tif",height=3,width=3,units = "in",res = 600,compression = "lzw")
s=c(22:31)+195
parameters=c(2:78)
counts=CP(s,parameters,min)
cornum$D2P3=counts
dev.off()
## 2nd phase T3 in D2
tiff("D2_T3.tif",height=3,width=3,units = "in",res = 600,compression = "lzw")
s=c(35:44)+195
parameters=c(2:78)
counts=CP(s,parameters,min)
cornum$D2P4=counts
dev.off()
## 2nd phase T4 in D2
tiff("D2_T4.tif",height=3,width=3,units = "in",res = 600,compression = "lzw")
s=c(49:58)+195
parameters=c(2:78)
counts=CP(s,parameters,min)
cornum$D2P5=counts
dev.off()
## 3rd phase in D2
tiff("D2_last.tif",height=3,width=3,units = "in",res = 600,compression = "lzw")
s=c(56:65)+195
parameters=c(2:78)
counts=CP(s,parameters,min)
cornum$D2P6=counts
dev.off()
#----------------------------------------
## 1st phase in D3
tiff("D3_1st.tif",height=3,width=3,units = "in",res = 600,compression = "lzw")
s=c(1:10)+260
parameters=c(2:78)
counts=CP(s,parameters,min)
cornum$D3P1=counts
dev.off()
## 2nd phase T1 in D3
tiff("D3_T1.tif",height=3,width=3,units = "in",res = 600,compression = "lzw")
s=c(6:15)+260
parameters=c(2:78)
counts=CP(s,parameters,min)
cornum$D3P2=counts
dev.off()
## 2nd phase T2 in D3
tiff("D3_T2.tif",height=3,width=3,units = "in",res = 600,compression = "lzw")
s=c(22:31)+260
parameters=c(2:78)
counts=CP(s,parameters,min)
cornum$D3P3=counts
dev.off()
## 2nd phase T3 in D3
tiff("D3_T3.tif",height=3,width=3,units = "in",res = 600,compression = "lzw")
s=c(35:44)+260
parameters=c(2:78)
counts=CP(s,parameters,min)
cornum$D3P4=counts
dev.off()
## 2nd phase T4 in D3
tiff("D3_T4.tif",height=3,width=3,units = "in",res = 600,compression = "lzw")
s=c(49:58)+260
parameters=c(2:78)
counts=CP(s,parameters,min)
cornum$D3P5=counts
dev.off()
## 3rd phase in D3
tiff("D3_last.tif",height=3,width=3,units = "in",res = 600,compression = "lzw")
s=c(56:65)+260
parameters=c(2:78)
counts=CP(s,parameters,min)
cornum$D3P6=counts
dev.off()
#--------------------------------------------------------------------------
View(cornum)
# output the numbers significant correlations per event per reactor
#--------------------------------------------------------------------------
#2.4.4 fate of subcommunities
abd=rbind(abundance[1,],abundance[2:66,],abundance[1,],abundance[67:131,],abundance[1,],abundance[132:196,],abundance[1,],abundance[197:261,],abundance[1,],abundance[262:326,])
abd$Group=c(rep("C1",66),rep("C2",66),rep("D1",66),rep("D2",66),rep("D3",66))
abd=abd[,2:71]
melt=melt(abd,id.vars = c("Group","Time_d"))

a=ggplot(data=melt[1:5610,],aes(x=Time_d,y=value))+facet_grid(variable~Group)+geom_line()+scale_x_continuous(name="Time(d)",breaks=c(0,45,90))+scale_y_continuous(name="Relative cells abundance per subcommunity",breaks=c(0,0.5),limits=c(0,0.5))+theme(axis.title = element_text(face="bold"))
a
b=ggplot(data=melt[5611:11220,],aes(x=Time_d,y=value))+facet_grid(variable~Group)+geom_line()+scale_x_continuous(name="Time(d)",breaks=c(0,45,90))+scale_y_continuous(name="Relative cells abundance per subcommunity",breaks=c(0,0.5),limits=c(0,0.5))+theme(axis.title = element_text(face="bold"))
b
c=ggplot(data=melt[11221:16830,],aes(x=Time_d,y=value))+facet_grid(variable~Group)+geom_line()+scale_x_continuous(name="Time(d)",breaks=c(0,45,90))+scale_y_continuous(name="Relative cells abundance per subcommunity",breaks=c(0,0.5),limits=c(0,0.5))+theme(axis.title = element_text(face="bold"))
c
d=ggplot(data=melt[16831:22440,],aes(x=Time_d,y=value))+facet_grid(variable~Group)+geom_line()+scale_x_continuous(name="Time(d)",breaks=c(0,45,90))+scale_y_continuous(name="Relative cells abundance per subcommunity",breaks=c(0,0.5),limits=c(0,0.5))+theme(axis.title = element_text(face="bold"))
d

#result show in a pdf file which named "FateofSubcommunities.pdf"
pdf("FateofSubcommunities.pdf",width=7.5,height = 10)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,2)))
print(a, vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(b, vp=viewport(layout.pos.row=1,layout.pos.col=2))

grid.newpage()
pushViewport(viewport(layout=grid.layout(1,2)))
print(c, vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(d, vp=viewport(layout.pos.row=1,layout.pos.col=2))
dev.off()
#===========================================================================
#End
#===========================================================================
