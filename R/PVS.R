print.PVS<-function(x,...){
cat("*** MLM fitted using PVS algorithm ***\n")
cat("\n-Call: ")
print(x$call)
cat("\n-Intercept (beta_0): ")
cat(x$beta.est[1])
cat("\n-Linear coefficients (beta): ")
cat(x$beta.est[-1])
cat("\n-Number of non-zero beta-coefficients: ")
cat(length(x$indexes.beta.nonnull[[1]]))
cat("\n-Indexes non-zero beta-coefficients: ")
cat(x$indexes.beta.nonnull[[1]])
cat("\n-Lambda: ")
cat(x$lambda.opt)
cat("\n-IC: ")
cat(x$IC)
cat("\n-Penalty: ")
cat(x$penalty)
cat("\n-Criterion: ")
cat(x$criterion)
cat("\n")
}


summary.PVS<-function(object,...){
x<-object
cat("*** MLM fitted using PVS algorithm ***\n")
cat("\n-Call: ")
print(x$call)
cat("\n-Intercept (beta_0): ")
cat(x$beta.est[1])
cat("\n-Linear coefficients (beta): ")
cat(x$beta.est[-1])
cat("\n-Number of non-zero beta-coefficients: ")
cat(length(x$indexes.beta.nonnull[[1]]))
cat("\n-Indexes non-zero beta-coefficients: ")
cat(x$indexes.beta.nonnull[[1]])
cat("\n-Lambda: ")
cat(x$lambda.opt)
cat("\n-IC: ")
cat(x$IC)
cat("\n-Penalty: ")
cat(x$penalty)
cat("\n-Criterion: ")
cat(x$criterion)
cat("\n")
}


predict.PVS<- function(object,newdata=NULL,y.test=NULL, ...)
{
if(is.null(newdata)){
	y <- fitted(object)
	out<-y
}
else{
	z.test <- newdata
	onez.test<- cbind(rep(1, nrow(z.test)), z.test)
	pred.LR <- as.matrix(onez.test)%*%object$beta.est
	MSEP<-  mean((pred.LR - y.test)^2)		
	if(is.null(y.test)) out<-pred.LR
	else out<-list(y=pred.LR,MSEP=MSEP)
}
out
}


plot.PVS<-function(x,ind=1:10,size=15,col1=1,col2=2,col3=4,option=0,...)
{
a<-x$range.grid[1]
b<-x$range.grid[2]
long_data <- as.data.frame(t(x$z[ind, ]))
Wavelength = seq(a, b, length.out = nrow(long_data))
long_data <- cbind(Wavelength, long_data)
long_data <- gather(long_data, key = "Series", value = "Value", -Wavelength)
matplot_colors <- 1:10
impact_points <- Wavelength[x$indexes.beta.nonnull[[1]]]
g1=ggplot(long_data, aes(x = Wavelength, y = Value, group = Series, color = Series)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = rep(matplot_colors, length.out = length(unique(long_data$Series)))) +
  geom_vline(xintercept = impact_points, color = 1, linetype = "dashed",linewidth=1) +  
  theme_bw() +
  labs(x = "t", y = "", title = expression(paste(zeta[i],"(t)"))) +
  theme(plot.title = element_text(hjust = 0.5,size=size),axis.title.x=element_text(size=size),legend.position = "none")  
fit<-fitted(x)
mod <- lm(x$y ~ fit)
res <- residuals(mod)
y<-x$y
x_df <- data.frame(fit, y)
g2<-ggplot(x_df, aes(x = fit, y = y)) +
  geom_point(colour = col1,shape=1, size=5) + 
  geom_smooth(method=lm,formula=y~x,colour =col2,  linewidth = 1.5) +
  labs(x = "Fitted values", y = "y", title = "Response vs Fitted values") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=size),axis.title.x=element_text(size=size))
x_df2 <- data.frame(fit, res)
g3<-ggplot(x_df2, aes(x = fit, y = res)) +
  geom_point(colour = col1,shape=1,size=5) + 
  geom_hline(yintercept = 0, linetype = "dashed", colour = 1, linewidth = 1) +
  geom_smooth(method=loess,formula=y~x,linewidth=1.5,col=col3)+
  labs(x = "Fitted values", y = "Residuals", title = "Residuals vs Fitted Values") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=size),axis.title.x=element_text(size=size))
if(option==0){
grid.arrange(g1,g2, g3,ncol =3)}
if(option==1){
grid.arrange(g1,ncol=1)}
if(option==2){
grid.arrange(g2,g3,ncol=2)}
}


