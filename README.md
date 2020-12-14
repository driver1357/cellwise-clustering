# cellwise-clustering
####Required libraries###########
# 

library("reshape2") # Para el melt
library("ggplot2")  # Para el ggplot
library("mvtnorm")
library("tclust")
library("splines")
library("robustbase")

#################################
#################################
####### Simulate data #########
#################################
#################################

set.seed(2)

#############################################################
######Simulation of data: scenario I(Spread contamination)  #########
#############################################################
# Para generaci?n de cellwise outliersa
alp.0 <- 0.03
# numero de grupos
K <- 2 
# Dimension
p<- 100
# Dimension intrinseca
q <- 2

# tama?os de los grupos (igual en los dos grupos)
n.i <- rep(200,K)
n <- sum(n.i[1:K])

# Generando datos (trabajamos en intervalo [0,1])
s <- (1:p)/p


m.1 <- function(x){5+10*sin(4*pi*x)*exp(-2*x)+5*sin(pi*x/3)+2*cos(2*pi/2)}
rho1.1 <- function(x){sqrt(2)*cos(2*pi*x)}
rho2.1 <- function(x){sqrt(2)*sin(2*pi*x)}

m.2 <- function(x){10+10*cos(4*pi*x)}
rho1.2 <- function(x){sqrt(2)*sin(2*pi*x)}
rho2.2 <- function(x){sqrt(2)*cos(2*pi*x)}

X <- matrix(NA,n,p)
datos.sin.error <- X
for (i in 1:n.i[1]){
    datos.sin.error[i,] <- m.1(s)+rnorm(1,0,3)*rho1.1(s)+rnorm(1,0,2)*rho2.1(s)
    var.error <- 0.5 
    X[i,] <- datos.sin.error[i,]+rnorm(p,0,var.error)
    
  }
  
for (i in (n.i[1]+1):(n.i[1]+n.i[2])){
    datos.sin.error[i,] <- m.2(s)+rnorm(1,0,4)*rho1.2(s)+rnorm(1,0,2)*rho2.2(s)
    var.error <- 0.5 
    X[i,] <- datos.sin.error[i,]+rnorm(p,0,var.error)
  }
  
  # Generando cellwise outliers dispersos (mitad y mitad a cada lado)
  so <- sample(n*p,n*p*alp.0)
  X[so] <- runif(n*p*alp.0,-20,-15)
  X[sample(so,floor(n*p*alp.0/2),replace=F)] <- runif(floor(n*p*alp.0/2),35,40)
  
  # R es la asignacion a clusters buena de las filas
  R<- rep(1:K, n.i[1:K])
  
  # Asignacion buena en Xk con 0 para celdas recortadas
  Xk <- matrix(NA,n,p)
  for (k in (1:K)){
    Xk[which(R==k),] <- k
  }
  Xk[so]<-0
  
  
#matplot(t(Ej.1[[2]]),type="l",lty=1)
matplot(t(X),type="l",lty=1)


#####################################################
#####################################################
######Algoritmo                             #########
#####################################################
#####################################################

# N. grupos
K <- 2
# Tamaño recorte
alp <- 0.05
# Dimension intrinseca
q <- 2

# Recorte trimmed k-means
alp.tk <- 0.1
# Recorte regresión robusta
alp.re <- 0.4
# Suavizado en lowess
f <- 1/5
# Número de nodos de splines
knots <- 4
# El número de loops 
loops <- 20

# Recorte final
alp.real <- alp.real



##################################
########Dimensionar matrices #####
x <- X
n<-nrow(X)
p<-ncol(X)
W<-array(NA,c(n,p,K))
W.reduc <-W
m <- matrix(0,nrow=p, ncol=K)
cls <- rep(1,n)
rij<-matrix(0,nrow=length(X), ncol=K)
rijk<-array(NA,c(n,p,K))
m<-matrix(0,nrow=p, ncol=K)
r<-matrix(0,nrow=n, ncol=K)
A<-array(NA,dim=c(q,n,K))
b<-array(NA,dim=c(q,p,K))
xhat<-array(NA,dim=c(n,p,K))

#####
# Generando datos (trabajamos en intervalo [0,1])
s <- (1:p)/p

##############################################################
#### PREPOCESADO ############################################
##############################################################

##################################
### Se suaviza robusto con "lowess"
lowess.mod <- function(x){lowess(x,f=f,iter=20)$y}
X.smooth <- t(apply(X,1,lowess.mod))

# Gráfico resumen curvas suavizadas
par(mfrow=c(1,1))
plot(s,X.smooth[1,],type="n",ylim=c(range(X)[1],range(X)[2]),xlab="",ylab="",main="Smoothed curves (lowess)")
for(i in 1:n){
  lines(s,X.smooth[i,],col=1)
}

##################################
### Se usa un ajuste de B-splines (suavizar más y reducir dimensionalidad para cluster
###   de los datos X.smooth)

# Número de nodos
knots <- knots
df <- knots+4

# XX será la representación en B splines
XX <- matrix(rep(0,n*df),nrow=n)
bb <- bs(s, df=df,intercept=T)
for (u in 1:n){
  fit <-lm(X.smooth[u,] ~ -1+bb)
  XX[u,] <-fit$coefficients
}

# Grafico resumen B.splines
par(mfrow=c(1,1))
plot(s,X[1,],type="n",ylim=c(range(X)[1],range(X)[2]),xlab="",ylab="",main="B-splines summary")
for(i in 1:n){
  lines(s,X[i,],col="grey")
}
for (i in 1:n){
  lines(s,bb[,]%*%XX[i,],col=i+1,pch=2)
}

##################################
# Trimmed k-means sobre XX

library(tclust)
for( i in 1:3){
  a <- tclust(XX,k=K,alpha=alp.tk,nstar=100*K,restr.fact=100) # Si se aumenta K hay que incrementar niter
  
}

# Gráfico del cluster
pairs(XX,col=a$cluster+1)

##################################
# Inicializar medias y autofunciones

## Medias
m <- bb[,]%*%a$centers[,]

## Autofunciones
for (k in 1:K){
  for (j in 1:q){
    b[j,,k] <- bb[,]%*%eigen(a$cov[,,k])$vectors[,j]
  }
}

## Graficos resumen inicialización cluster
par(mfrow=c(1,1))
plot(s,X[1,],type="n",ylim=c(range(X)[1],range(X)[2]),xlab="",ylab="",main="Initial clusters")
for(i in 1:n){
  lines(s,X[i,],col="grey")
}
for (k in 1:K){
  ind <- which(a$cluster==k)
  for (i in ind){
    lines(s,bb[,]%*%XX[i,],col=k+1)
  }
}
for (k in 1:K){
  lines(s,m[,k],lwd=2,lty=2)
}

# Graficos resumen inicialización autofunciones
par(mfrow=c(k,q))
for (k in 1:K){
  for (j in 1:q){
    plot(s,b[j,,k],xlab="",ylab="",main=paste("k=",k,"and","q=",j),ylim=c(range(b)[1],range(b)[2])) 
    abline(h=0,col=2)
  }
}

##################################
### Obtener primeros A's, W y W.reduc con regresión robusta y primera asignación
###  (esta parte es lenta: se podría hacer regresión normal y definir r[i,k] con la suma de los 60% residuales más pequeños)

for (k in (1:K)){
  #Datos centrados
  xc<-scale(x,center=m[,k],scale = FALSE )
  
  #Estimacion de A (era un poco más lioso antes porque la respuesta es identicamente nula media igual a una observacion)
  for(i in 1:n) {
    rr <- ltsReg(xc[i,]~matrix(t(b[1:q,1:p,k]),nrow=p,ncol=q)-1,alpha=1-alp.re)   # 1 - 0.6 = 0.4 es un recorte grande
    A[,i,k] <-rr$coefficients
    r[i,k] <- rr$raw.scale  #Usa los residuales para crear una distancia a clusters
  }
}

# Clasificación resultante del LTS
cls<-apply(r,1,which.min)

#X estimada
for(k in 1:K){
  for(i in 1:n){
    xhat[i,,k]<-m[,k]+as.vector(matrix(t(b[1:q,1:p,k]),nrow=p,ncol=q)%*%(matrix(A[1:q,i,k],nrow=q,ncol=1)))
    rijk[i,,k]<-(x[i,]-xhat[i,,k])^2
  }
}

# Crear W
n.cls <- table(cls)
for(k in 1:K){
  for (k1 in 1:K){
    for(j in 1:p){
      W[cls==k1,j,k] <- as.numeric( rijk[cls==k1,j,k] <= sort(rijk[cls==k1,j,k])[ceiling(n.cls[k1]*(1-alp))] )
    }
  }
  r[,k]<-(p/apply(W[,,k],1,sum))*apply(W[,,k]*rijk[,,k],1,sum)
}

# Nueva asignación con W  (esta parte se puede no hacer y dejar fijo algún "cls" anterior)
cls<-apply(r,1,which.min)
print(table(cls))

#  Crear W.reduc
for(k in 1:K){
  W.reduc[,,k] <- W[,,k]
  W.reduc[which(cls!=k),,k]<-rep(0,p)
}

# heatmap(rijk[,,2],scale="none",Colv=NA,Rowv=NA)
# heatmap(W[,,2],scale="none",Colv=NA,Rowv=NA)
# heatmap(W.reduc[,,2],scale="none",Colv=NA,Rowv=NA)


##################################
### Iteraciones (loops) ##########

# Inicializar shat.opt
shat.opt <- 10^10

# Numero de loops
loops <- loops  

#inicia el loop
it<-0

# Todo va bien si abortar es 1
abortar <- 1

while((it<loops)&(abortar==1)){
  it<-it+1
  print(it)
  
  for(k in 1:K){
    ###Estimacion en b
    rida<-apply((W.reduc[,,k]),2,sum)
    ida<-which(rida>(q+1))
    for(j in sort(ida)){
      b[1:q,j,k] <- lm((x[,j]-m[j,k]) ~matrix(t(A[1:q,1:n,k]),nrow=n,ncol=q)-1,weights = W.reduc[,j,k])$coefficients
    }
  }
  
  ### estimacion de la media
  rm <- X
  for(k in 1:K){
    
    for(i in 1:n){
      rm[i,1:p]<-x[i,]-as.vector(matrix(t(b[1:q,1:p,k]),nrow=p,ncol=q)%*%(matrix(A[1:q,i,k],nrow=q,ncol=1)))
    }
    
    m[,k]<-as.matrix((1/colSums(W.reduc[,,k]))*as.matrix(colSums(W.reduc[,,k]*rm)))
  }
  
  ######
  for(k in 1:K){
    # Datos centrados
    xc<-scale(x,center = m[,k],scale = FALSE)
    
    ###Estimacion en A
    iw<-apply(as.matrix(W[,,k]),1,sum)
    inx<-which(iw>(q+1))
    for(i in inx) {
      A[1:q,i,k] <- lm(xc[i, ]~matrix(t(b[1:q,1:p,k]),nrow=p,ncol=q)-1,weights = W[i,1:p,k])$coefficients
    }
  }
  
  ###X estimada
  for(k in 1:K){
    for(i in 1:n){
      xhat[i,,k]<-m[,k]+as.vector(matrix(t(b[1:q,1:p,k]),nrow=p,ncol=q)%*%(matrix(A[1:q,i,k],nrow=q,ncol=1)))
      rijk[i,,k]<-(x[i,]-xhat[i,,k])^2
    }
  }
  
  ## W con todos
  n.cls <- table(cls)
  for(k in 1:K){
    for (k1 in 1:K){
      for(j in 1:p){
        W[cls==k1,j,k] <- as.numeric( rijk[cls==k1,j,k] <= sort(rijk[cls==k1,j,k])[ceiling(n.cls[k1]*(1-alp))] )
      }
    }
    r[,k]<-(p/apply(W[,,k],1,sum))*apply(W[,,k]*rijk[,,k],1,sum)
  }
  
  
  # Assignments  (se podría dejar fijo y usar "cls" inicialización comentando la siguiente linea)
  cls<-apply(r,1,which.min)
  print(table(cls))
  
  ## W.reduc
  for(k in 1:K){
    W.reduc[,,k] <- W[,,k]
    W.reduc[which(cls!=k),,k]<-rep(0,p)
  }
  
  ## Ver si todo va bien
  chequeo <- rep(0,K)
  for (k in (1:K)){
    chequeo[k]=sum(W[which(cls==k),,k])
  }
  #print(chequeo)
  
  # Calcula la funcion objetivo o aborta si tiene problemas con grupos vacios
  if (min(chequeo)>10){shat <- 0
  for(k in 1:K){
    if( is.na( sum(rijk[,,k]) ) ==F )   {shat <- shat +sum(W.reduc[,,k]*rijk[,,k])}  else {shat <- shat +10^10}
  }
  print(shat)
  }
  else{abortar <-0
  print("se aborto")
  #print(chequeo)
  shat <- 10^10
  }
  
}  # Fin del loop


##################################
### Resultados          ##########

W.opt <- W
cls.opt <- cls
print(table(cls,rep(1:K,n.i[1:K])))  # Dara error si no hay n.i
shat.opt <- shat
A.opt <- A
b.opt <- b
print(shat.opt)
m.opt <- m
rijk.opt <- rijk

# Preciciones en array nxpxK
xhat.opt <- xhat

# Asisgnaciones individuales en WC (matriz nxp con cada celda a su grupo y 0 si es recortada)
WC<-matrix(0,n,p)
ig<-which(cls.opt==1)
WC[ig,]<-(W.opt[ig,,1])
for(i in 2:K){
  ig<-which(cls.opt==i)
  Wk<-W.opt[,,i]
  id<-which(Wk==1)
  Wk[id]<-i
  WC[ig,]<-Wk[ig,]
}

### Predicciones en matriz nxp
xpred <- x
for (i in 1:n){
  xpred[i,]<-xhat.opt[i,,cls.opt[i]]
}


##################################
##################################
### Gráficos resultados ##########
##################################

# Grafico de lineas grupos juntos y circulos en recortadas
par(mfrow=c(1,1))
plot(s,X[1,],type="n",col=cls.opt[1]+1,ylim=c(range(X)[1],range(X)[2]),xlab="",ylab="")
for(i in 1:n){
  lines(s,X[i,],col=cls.opt[i]+1)
}

for(i in 1:n){
  points(s[which(W.opt[i,,cls.opt[i]]==0)],X[i,W.opt[i,,cls.opt[i]]==0],cex=0.6)
}


# Grafico de lineas por grupos separados y circulos en recortadas
par(mfrow=c(2,ceiling(K/2)))
for (k in (1:K)){
  X.grupo <- X[cls.opt==k,]
  W.grupo <- W.opt[cls.opt==k,,k]
  plot(s,X.grupo[1,],type="n",col=k+1,ylim=c(range(X)[1],range(X)[2]),xlab="",ylab="",main=paste("Cluster",k))
  for (i in (1:dim(X.grupo)[1])){
    lines(s,X.grupo[i,],col=k+1)
  }
  for (i in (1:dim(X.grupo)[1])){
    points(s[which(W.grupo[i,]==0)],X.grupo[i,W.grupo[i,]==0],cex=0.6)
  }
  lines(s,m.opt[,k],lwd=3)
}

# # Grafico grupos separados con las coordenadas no recortadas
# par(mfrow=c(1,K))
# for (k in (1:K)){
#   X.grupo <- X[cls.opt==k,]
#   W.grupo <- W.opt[cls.opt==k,,k]
#   plot(s,X.grupo[1,],type="n",col=k+1,ylim=c(range(X)[1],range(X)[2]),xlab="",ylab="",main=paste("Cluster",k))
#   for (i in (1:dim(X.grupo)[1])){
#     points(s[which(W.grupo[i,]!=0)],X.grupo[i,W.grupo[i,]!=0],pch=19,cex=0.5)
#   }
#   lines(s,m.opt[,k],col=k+1,lwd=3)
# }

# Grafico autofunciones estimadas
par(mfrow=c(K,q))
for (k in 1:K){
  for (j in 1:q){
    plot(s,b.opt[j,,k],xlab="",ylab="",main=paste("k=",k,"and","q=",j),ylim=c(range(b)[1],range(b)[2]),type="b") 
    abline(h=0,col=2)
  }
}

# Grafico para interpretar variabilidades
par(mfrow=c(K,q+1))
for (k in 1:K){
  X.grupo <- X[cls.opt==k,]
  W.grupo <- W.opt[cls.opt==k,,k]
  plot(s,X.grupo[1,],type="n",col=k+1,ylim=c(range(X)[1],range(X)[2]),xlab="",ylab="",main=paste("Cluster",k))
  for (i in (1:dim(X.grupo)[1])){
    lines(s,X.grupo[i,],col=k+1)
  }
  
  for (j in 1:q){
    const <- sqrt(var(A.opt[j,cls.opt==k,k]))
    plot(s,m.opt[,k]+b.opt[j,,k],xlab="",ylab="",main=paste("k=",k,"and","q=",j),ylim=c(range(X)[1],range(X)[2]),type="n") 
    points(s,m.opt[,k]+2*const*b.opt[j,,k],xlab="",ylab="",pch="+",ylim=c(range(X)[1],range(X)[2]),col=2)
    points(s,m.opt[,k]-2*const*b.opt[j,,k],xlab="",ylab="",pch="-",ylim=c(range(X)[1],range(X)[2]))
    
  }
}

# Grafico de los "scores"
if (q>=2){
  par(mfrow=c(2,ceiling(K/2)))
  for(k in 1:K){
    # AA <- A.opt[1:2,cls.opt==k,k]
    # plot(t(A.opt[1:2,cls.opt==k,k]),main=paste("Scores",k),xlab="",ylab="",type="n",xlim=c(range(AA)[1],range(AA)[2]),ylim=c(range(AA)[1],range(AA)[2]))
    plot(t(A.opt[1:2,cls.opt==k,k]),main=paste("First two scores of cluster ",k),xlab="",ylab="",type="n")
    text(t(A.opt[1:2,cls.opt==k,k]),cex=0.8)
  }
}

# Gráfico bonito de WC (se pueden poner nombres a individuos)
colnames(WC) <- seq(1,p)
rownames(WC) <- nombres
longData<-melt(WC)
ggplot(longData, aes(x = Var2, y = Var1, fill=as.factor(value),color=as.factor(value))) +
  geom_raster() +   labs(x="p", y="n")+
  scale_fill_manual(name="code",values=1:(K+1) )




### Comparación predicciones y observado
# Solo las no-recortadas
heatmap((x-xpred)*(WC==0),scale="none",Colv=NA,Rowv=NA,col=rainbow(256,start = 1/6))
# Todas
heatmap((x-xpred),scale="none",Colv=NA,Rowv=NA,col=rainbow(256,start = 1/6))
heatmap((x-xpred),scale="none",Colv=NA,Rowv=NA,col=cm.colors(256))

### Más atipicas
mas <- rep(0,n)
for (i in 1:n){
  mas[i]<-sum((x-xpred)[i,]^2*(WC[i,]==0))
}
# Cuales?
sort(mas,decreasing = TRUE,index=TRUE)$ix
nombres[sort(mas,decreasing = TRUE,index=TRUE)$ix[1:10]]

# Grafico
par(mfrow=c(1,1))
plot(1:n,mas,type="b",xlab="index",ylab="sum of squared residuals")

#### Algún grafico individual (discontinuo es la predicción)
ii <- 3
par(mfrow=c(1,1))
plot(s,x[1,],type="n",ylim=c(range(X)[1],range(X)[2]),xlab="",ylab="",main=paste(nombres[ii]))
for (i in 1:n){
  lines(s,x[i,],col="grey")
}
lines(s,x[ii,])
lines(s,xpred[ii,],lty=2,lwd=3)
points(s[which(WC[ii,]==0)],X[ii,WC[ii,]==0],cex=0.6)

#### Algunos grafico individuales (6 al azar)
par(mfrow=c(2,3))
iii <- sample(1:n,6,replace=TRUE)
for (ii in iii){
  plot(s,x[1,],type="n",ylim=c(range(X)[1],range(X)[2]),xlab="",ylab="",main=paste(nombres[ii])) # Se podria usar el nombre del individuo
  for (i in 1:n){
    lines(s,x[i,],col="grey")
  }
  lines(s,x[ii,])
  lines(s,xpred[ii,],lty=2,lwd=3)
  points(s[which(WC[ii,]==0)],X[ii,WC[ii,]==0],cex=0.6)
}

#### Algunos grafico individuales (6 más atípicos)
par(mfrow=c(2,3))
iii <- sort(mas,decreasing = TRUE,index=TRUE)$ix[1:6]
for (ii in iii){
  plot(s,x[1,],type="n",ylim=c(range(X)[1],range(X)[2]),xlab="",ylab="",main=paste(nombres[ii])) # Se podria usar el nombre del individuo
  for (i in 1:n){
    lines(s,x[i,],col="grey")
  }
  lines(s,x[ii,])
  lines(s,xpred[ii,],lty=2,lwd=3)
  points(s[which(WC[ii,]==0)],X[ii,WC[ii,]==0],cex=0.6)
}


#############################################################
###### Posibilidad de refininamiento                #########
#############################################################

# Los residuales al cuadrado de las observaciones
Rij<-matrix(0,ncol=p, nrow=n)
for(i in 1:length(cls.opt)){
  Rij[i,]<-rijk.opt[i,1:p,cls.opt[i]]
}

# Dar alp.real si lo conocemos
alp.real <- alp.real
# Gráfico para estimar el tamaño de recorte final con los nxp residuales ordenados decrecientemente (linea vertical es el alp.real)
par(mfrow=c(1,1))
plot((1:ceiling(n*p*0.1))/(n*p),sort(Rij,decreasing = TRUE)[1:ceiling(n*p*0.1)],ylab="residual",xlab="prop")
abline(v=alp.real)

# Si conocieramos el alp.real (que se podría hacer estudiando la grafica anterior...)
#   podemos sacar WC.final con la asignación final
Rijo.out<-sort(Rij,index.return=TRUE,decreasing = TRUE)$ix[1:(alp.real*n*p)]
WC.final<-matrix(0,nrow=n, ncol=p)
for(k in 1:K){
  WC.final[which(cls==k),]<-rep(k,p)
}
WC.final[Rijo.out] <- 0

## Curvas recortadas enteras
ttt <- which(apply(WC.final==0,1,sum)/p > alp.re)
nombres[ttt]

### Cambiar WC.final quitando recortadas enteras
for(i in ttt){
  WC.final[i,] <- rep(0,p)
}

# Grafico bonito final de WC.final
colnames(WC.final) <- seq(1,p)
rownames(WC.final) <- nombres  # O el nombre de los individuos nombres
longData.final<-melt(WC.final)
ggplot(longData.final, aes(x = Var2, y = Var1, fill=as.factor(value),color=as.factor(value))) +
  geom_raster() +   labs(x="p", y="n")+
  scale_fill_manual(name="code",values=1:(K+1) )

### Comparación predicciones y observado
# Solo las no-recortadas y muy roja la recortada enteramente
represent <- (x-xpred)*(WC.final==0)
maxx <- range(represent)[2]
for (i in ttt){
  represent[i,] <- rep(maxx,p)
}
heatmap(represent,scale="none",Colv=NA,Rowv=NA,col=rainbow(256,start = 1/6))

### Comparación predicciones y observado
# y muy roja la recortada enteramente
represent <- (x-xpred)
maxx <- range(represent)[2]
for (i in ttt){
  represent[i,] <- rep(maxx,p)
}
heatmap(represent,scale="none",Colv=NA,Rowv=NA,col=rainbow(256,start = 1/6))

# Gráfico de medias
par(mfrow=c(1,1))
plot(s,X[1,],type="n",col=cls.opt[1]+1,ylim=c(range(X)[1],range(X)[2]),xlab="",ylab="")
for(i in 1:n){
  lines(s,X[i,],col="grey")
}
for (k in (1:K)){
  lines(s,m.opt[,k],lwd=2,col=k+1)
}

# Grafico por grupos por seaparados despues refinamiento
par(mfrow=c(2,ceiling(K/2)))
for (k in (1:K)){
  plot(s,X[1,],type="n",col=k+1,ylim=c(range(X)[1],range(X)[2]),xlab="",ylab="",main=paste("cluster",k))
  cluster <- which(cls.opt==k)
  for (i in cluster){
    lines(s,X[i,],col=k+1)
  }
  #lines(s,m.opt[,k],lwd=3)
  for (i in setdiff(cluster,ttt)){
    points(s[which(WC.final[i,]==0)],X[i,WC.final[i,]==0],cex=0.7)
  }
  if (length(intersect(cluster,ttt))>0){
    for (ii in intersect(cluster,ttt)){
      points(s,X[ii,],lwd=2,cex=0.7)
    }
  }
}


# Trimmed only
par(mfrow=c(1,1))
plot(s,X[1,],type="n",col=cls.opt[1]+1,ylim=c(range(X)[1],range(X)[2]),xlab="",ylab="")
for(i in 1:n){
  lines(s,X[i,],col="gray")
}

for(i in ttt){
  lines(s,X[i,],col=1,cex=0.7,lwd=2)
}
for(i in setdiff(1:n,ttt)){
  points(s[which(WC.final[i,]==0)],X[i,WC.final[i,]==0],cex=0.7)
}




# Grafico de lineas grupos juntos y recortadas despues refinamiento
par(mfrow=c(1,1))
plot(s,X[1,],type="n",col=cls.opt[1]+1,ylim=c(range(X)[1],range(X)[2]),xlab="",ylab="")
for(i in 1:n){
  lines(s,X[i,],col=cls.opt[i]+1)
}
for(i in ttt){
  points(s,X[i,],col=1,lwd=3,cex=0.7)
}
for(i in setdiff(1:n,ttt)){
  points(s[which(WC.final[i,]==0)],X[i,WC.final[i,]==0],cex=0.7)
}

