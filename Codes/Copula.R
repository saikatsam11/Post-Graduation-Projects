library(fitdistrplus)
dd=data(danishmulti)
dd

str(danishmulti)
x_data=data.frame(danishmulti$Building,danishmulti$Contents)
x_data[,1]=danishmulti$Building
cor(x_data$danishmulti.Building,x_data$danishmulti.Contents)
dim(danishmulti$Building)
x_data[,1]=as.numeric(x_data[,1])
summary(x_data[,1])

x_data$danishmulti.Building
class(danishmulti$Building)

write.csv(danishmulti,"saikatdanish")
get.wd()
getwd()
x_data[1]=as.numeric(x_data[1])
x_data[2]=as.numeric(x_data[2])

######################################################################
data12=read.csv(file.choose(),header = T)
str(data12)
cor(data12$Building,data12$Content)
#########################

x_data=data.frame(danishmulti$Building,danishmulti$Contents)

ind <- which(x_data$danishmulti.Building >0 & x_data$danishmulti.Contents >0  )
summary(ind)

claim.nozero <- x_data[ind,]
summary(claim.nozero)
cor(claim.nozero[,1],claim.nozero[,2],method = "kendall")
cor(claim.nozero[,1],claim.nozero[,2])
length(claim.nozero[,1])
quantile(claim.nozero$danishmulti.Building)
plot(log(claim.nozero[,1]),log(claim.nozero[,2]))
hist(claim.nozero[,1], breaks = 130, main = '', xlab = 'Loss to Buildings',xlim = c(0,95))
hist(claim.nozero[,2],,breaks = 130,main = '',xlab = 'Loss to Contents')
library(plot3D)
z <- table(claim.nozero[,1], claim.nozero[,2])
hist3D(z=z, border="black")

library(VGAM); 
library(stats4);            #for the function 'mle''

# --------------- defining exponential and exponentially tempered Pareto distribution ----------- #


pexponential = function(x, para){
  lmd = para[1];
  return( ifelse( x>0, (lmd*exp(-(lmd*x))), 0 ) );  
}

dexponential = function(x, para){
  lmd = para[1];
  return( ifelse( x>0, (1 - exp(-(lmd*x))) , 0 ) );   
}

pexptempareto = function(x,para){
  alpha = para[1];  beta = para[2];  theta = para[3];
  return( ifelse( x>0, ((theta^alpha) * exp(beta * theta) * x^(-(alpha + 1)) * exp(-(beta * x)) * (alpha + (beta * x))) , 0 ));
}

dexptempareto = function(x,para){
  alpha = para[1];  beta = para[2];  theta = para[3];
  return( ifelse( x>0, (1 - ((theta^alpha)*exp(beta * theta) * x^(-(alpha)) * exp(-(beta * x)))) , 0 ))
}




Quietmode = "ON";
#------------------------------------------------------
#Main computation functions
f1_truncated = function(x, para, splitpt){
  return ( f1(x, para)/Bigf1(splitpt, para) );
}
f2_truncated = function(x, para, splitpt){
  return ( f2(x, para)/(1-Bigf2(splitpt, para)) );
}
composite_f = function(x, para1, para2, r_theta){
  #The composite has parameters #para1 + #para2 + 2 (r and theta) in total
  #Continuity condition gives expression of r, so we have #para1 + #para2 + 1
  r = r_theta[1];
  splitpt = r_theta[2];
  if (x < splitpt){
    return ( r * f1_truncated(x, para1, splitpt) );
  }else{
    return ( (1-r) * f2_truncated(x, para2, splitpt) );
  }
}

composite_Bigf = function(x, para1, para2, r_theta){
  #The composite has parameters #para1 + #para2 + 2 (r and theta) in total
  #Continuity condition gives expression of r, so we have #para1 + #para2 + 1
  r = r_theta[1];
  splitpt = r_theta[2];
  if (x < splitpt){
    return ( r * Bigf1(x,para1)/Bigf1(splitpt, para1) );
  }else{
    return ( r + (1-r) * (Bigf2(x,para2)-Bigf2(splitpt,para2))/(1-Bigf2(splitpt, para2)) );
  }
}

ETP_complete_para = function(para){
  #Remapping
  alpha = para[1]; beta = para[2]; theta= para[3];
  #Complete the parameter set
  splitpt = theta;        
  if(splitpt<=0){
    lmd =NA;
  }else{
    lmd = (1/theta)*(alpha + (alpha + (theta*beta))^2);
  }  
  #Remapping  
  para1 = c(lmd);    para2 = c(alpha, beta, theta);
  #Compute the mixture parameters, the mixing weight and the splitpt.
  r = f2(splitpt, para2)*Bigf1(splitpt, para1) / (f2(splitpt, para2)*Bigf1(splitpt, para1) + f1(splitpt, para1)*(1-Bigf2(splitpt, para2)) );  
  r_splitpt = c(r, splitpt);
  #Return list of 3 vectors.  
  return( list(para1, para2, r_splitpt) );
}
ETP_para_condition = function(para1, para2, r_splitpt){  
  #Overall NA check.    
  if (is.all_not_na(c(para1,para2,r_splitpt))==FALSE){return(FALSE);}
  
  #Now value check.
  #condition for parameters in the first distribution
  lmd = para1[1];
  
  #condition for parameters in the second distribution
  alpha = para2[1];  beta = para2[2];  theta = para2[3];
  
  #condition for r and splitpt
  r = r_splitpt[1]; splitpt = r_splitpt[2];
}


name_list = c("pexptempareto","dexptempareto", "pexponential", "dexponential", "ETP_para_condition","ETP_complete_para");


#Assignment of functions
f1=getFunction(name_list[3]);      Bigf1=getFunction(name_list[4]);
f2=getFunction(name_list[1]);      Bigf2=getFunction(name_list[2]);
para_condition = getFunction(name_list[5]);
complete_para = getFunction(name_list[6]);


estimation=function(x_vec)
{
  NLL = function(para, x_vec){  
    cp = complete_para(para);  
    para1 = as.vector( cp[[1]] );
    para2 = as.vector( cp[[2]] );
    r_splitpt = as.vector( cp[[3]] );  
    
    #  if( !para_condition(para1, para2, r_splitpt) ){ return(NA) };
    
    x_vec = sort(x_vec);  
    L = 0;
    for (x in x_vec){           
      L = L + log( composite_f(x, para1, para2, r_splitpt));   
    }
    return(-L);
  }
  
  xx1=nlm(NLL,x_vec=claim.nozero[,1],c(1.733798699,0.005574184,7.975004616),hessian = TRUE) 
  xx2=nlm(NLL,x_vec=claim.nozero[,2],c(0.95672468 ,0.02348697,1.22266390),hessian = TRUE)
 
  parameter1= ETP_complete_para(para=c(xx1$estimate[1],xx1$estimate[2],xx1$estimate[3],xx1$estimate[4]))
  parameter2= ETP_complete_para(para=c(xx2$estimate[1],xx2$estimate[2],xx2$estimate[3],xx2$estimate[4]))
  para11=as.vector(parameter1[[1]])
  para12=as.vector(parameter1[[2]])
  r_theta1=as.vector(parameter1[[3]])
  
  para21=as.vector(parameter2[[1]]) 
  para22=as.vector(parameter2[[2]])
  r_theta2=as.vector(parameter2[[3]])
  
  par <- list(para11 = para11,para12 = para12, r_theta1 = r_theta1,para21=para21,para22=para22,r_theta2=r_theta2)
  
  n=length(claim.nozero[,1])
  # The fitted CDF
  #u1 <- sapply(1:n,function(k) composite_Bigf(claim.nozero[,1][k], para11, para12, r_theta1)) #MTPL
  #u2 <- sapply(1:n,function(k) composite_Bigf(claim.nozero[,2][k], para21, para22, r_theta2)) #MTPL
  
  u1 <- sapply(1:n,function(k) composite_Bigf(ss[,1][k], para11, para12, r_theta1))  ##DANISH
  u2 <- sapply(1:n,function(k) composite_Bigf(ss[,2][k], para21, para22, r_theta2))  ##DANISH
  return(c(par,list(u_fit =cbind(u1,u2),value = c(xx1$minimum,xx2$minimum))))
  
}

ss=cbind(claim.nozero[,1],claim.nozero[,2])
est=estimation(ss)
est

library(rvinecopulib)
library(cramer)

estimation_cp <- function(u_fit){
  fit_cop <- bicop(u_fit, family_set = 'joe',var_types = c("c", "c"))
  return(list(par = fit_cop$parameters,value = c(fit_cop$loglik)))
}


copula_parameter=estimation_cp(est$u_fit)
fit_cop <- bicop(est$u_fit, family_set = 'joe',var_types = c("c", "c"))
summary(fit_cop)
contour(fit_cop,col="red")

summary_model <- function(mar_fit,cop_fit,nobs){
  npar <- 7
  loglike <- sum(mar_fit$value) - cop_fit$value
  AIC <- 2*loglike + 2*npar
  BIC <- 2*loglike + log(nobs)*npar
  return(c(loglike,AIC,BIC))
}
summary_model(est,copula_parameter,nrow(claim.nozero))
sum(est$value)
sum(est$value) - copula_parameter$value

library(VineCopula)
BiCopPar2Tau(family = 6 , par = 1.4)
BiCopPar2TailDep(family = 6 , par = 1.4) 
