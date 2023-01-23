semimetric.projec<- function(data1, data2, 
theta,  order.Bspline=3, nknot.theta=3,
range.grid=NULL,nknot=NULL)
{
if(is.vector(data1)) data1 <- as.matrix(t(data1))
if(is.vector(data2)) data2 <- as.matrix(t(data2))
testfordim <- sum(dim(data1)==dim(data2))==2
twodatasets <- T
if(testfordim) twodatasets <- sum(data1==data2)!=prod(dim(data1))
projec1 <-projec(data=data1,theta=theta, range.grid=range.grid,order.Bspline=order.Bspline,nknot=nknot,nknot.theta=nknot.theta)
if(twodatasets){
	projec2 <-projec(data=data2,theta=theta, range.grid=range.grid,order.Bspline=order.Bspline,nknot=nknot,nknot.theta=nknot.theta)
}
else{
	projec2 <- projec1
}
SEMIMETRIC <- outer(projec1,projec2,"-")
table<-as.matrix(abs(SEMIMETRIC))
r<-nrow(table)
v<-ncol(table)
rownames(table)<-paste("D1",1:r,sep=".")
colnames(table)<-paste("D2",1:v,sep=".")
table
}

