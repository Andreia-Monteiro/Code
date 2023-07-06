#####2023-06-29
#Functions needed for Construct Covariate for Preferential Sampling paper
#Andreia Monteiro
#Isabel Natário



simul.rand.effect<-function(range0,sigma0,mesh,mu.field=4){
  #Function to simulate the spatial random effect with mean mu.field
  kappa0 <- sqrt(8)/range0     
  tau0 <- 1/(sqrt(4*pi)*kappa0*sigma0)
  spde_matern <- inla.spde2.matern(mesh=mesh, B.tau=cbind(log(tau0), -1, +1),
                                   B.kappa=cbind(log(kappa0), 0, -1), theta.prior.mean=c(0,0),
                                   theta.prior.prec=c(0.1,0.1))  #eventualmente podemos permitir que estes hyperparametros sejam fornecidos pelo utilizador
  
  Q <- inla.spde.precision(spde=spde_matern, theta=c(0,0))  #idem
  u <- inla.qsample(n=1, Q=Q, mu=rep(mu.field,nrow(Q))) #média mu.field;
  
  # Simulation grid
  
  limlattice <- c(0,1) # limits of the lattice simulation - unit square domain
  grid <- as.matrix(expand.grid(
    seq(limlattice[1], limlattice[2], length.out=100),
    seq(limlattice[1], limlattice[2], length.out=100)))
  
  # Simulated random field
  u_ver <- as.vector(inla.spde.make.A(mesh=mesh, loc=grid)%*%as.vector(u))
  
  return(list(u_ver=u_ver,grid=grid))
}


plot.simulated.field<-function(grid,u_ver){
#This function plots the simulated S(x)
  
ggplot() + geom_tile(data=data.frame(x=grid[,1], y=grid[,2],u=u_ver), 
                     mapping=aes(x=x,y=y,fill=u)) +
  scale_fill_viridis_c(option = "B",direction = -1) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  labs(fill="S(x)")
}


sample.prefer<-function(u_ver,nugget.prec=1/0.1,grid,n,beta=2){
  #This function simulate data preferentially according to a PP with intensity exp(beta*S(x))
  
  #Simulating a gaussian response variable
  predictor <- u_ver
  prec.gaussian <- nugget.prec  #nugget precision
  
  ysim <- rnorm(nrow(grid), mean=predictor, sd=prec.gaussian**(-1/2))
  DFSim <- data.frame(x=grid[,1], y=grid[,2], ysim=ysim)
  
  # Sample the data according to PP with intensity exp(beta*S(x))
  sampData <- DFSim[sample(1:nrow(DFSim), size=n, prob = exp(beta*u_ver)),] 
  
  return(sampData)

}


MLC<-function(domain_boundary,sampData){
  #MLC test data for preferentiability
  
  x.aux <- st_sf(st_sfc(st_polygon(list(domain_boundary))))
  
  #Test grid size
  dat<-sampData
  names(dat)<-c("x","y","resp")  
  
  #Number of points in sample
  n<-dim(dat)[1]
  
  #Final decision: consider the point density
  area<-st_area(x.aux)
  grid.size<-sqrt(area/n)  #size grid given by the square root of the area divided by the number of points (density)
  #grid.size might be a vector of different sizes, if defined otherwise; programming is made
  #considering that grid.size might be a vector of values

  res.1<-data.frame(grid.size=grid.size,corr=NA,p.value=NA)
  
  for (lg in 1:length(grid.size)){
    
    #grid
    grid <- st_make_grid(x.aux, cellsize = c(grid.size[lg],grid.size[lg]))
    
    #Determine which grid cells are inside the polygon: those with all the vertices grid[[i]]
    #inside domain_boundary
    #library(mgcv)
    grid.in<-NULL
    for (i in 1:length(grid)){
      aux<-in.out(as.matrix(domain_boundary),as.matrix(grid[[i]])) # true (=1) if any point of grid[[i]] is inside the domain_boundary
      if(sum(aux)==5) #5-1 grid cell vertices; last point is to close the cell
      {grid.in[i]<-"inside"}
      else{
        if(sum(aux)==0) grid.in[i]<-"fora" else grid.in[i]<-"border"}
    } #end for i
    
    #Data frame to keep the information for the test
    data.frame.test<-data.frame(ID=1:length(grid.in),grid.pos=grid.in,no.points=NA,y.average=NA)
    
    #Which grid cells that are totally inside the domain_boundary
    for (i in 1:length(grid)){
      if (data.frame.test$grid.pos[i]=="inside" | data.frame.test$grid.pos[i]=="border"){
        aux<-in.out(as.matrix(grid[[i]]),as.matrix(data.frame(dat$x,dat$y))) #Which data points are inside the cell grid
        data.frame.test$no.points[i]<-sum(aux)
        if (sum(aux)!=0) data.frame.test$y.average[i]<-mean(dat$resp[aux],na.rm=T)
        else data.frame.test$no.points[i]<-NA
        # cat(paste("Iteration:",i,"\n"))
      }
    }
    data.frame.test<-data.frame.test[data.frame.test$grid.pos=="inside"|data.frame.test$grid.pos=="border",]
    
    test.result<-cor.test(data.frame.test$no.points, data.frame.test$y.average,  method = "spearman", use="pairwise.complete.obs",exact=FALSE)
    res.1$corr[lg]<-test.result$estimate
    res.1$p.value[lg]<-test.result$p.value
  } #end for lg
  
  return(res.1)
}


construct.covariate<-function(sampData,perc=5){
#Function to obtain a constructed covariate given by the average distance to the k 
#nearest neighbours (percx100% nearest neighbours)

x<-sampData[,1]
y<-sampData[,2]
marks<-sampData[,3]
n<-dim(sampData)[1]

viz<-(perc/100)*n
nndistance1<-nndist(x,y, k=1:viz)  #computes the distances to the k=viz nearest neighbours
nndistance2<-as.data.frame(nndistance1)
nndistance<-apply(nndistance2,1,mean)

return(nndistance)
}


plot.covariate<-function(sampData,nndistance,grid,u_ver){
  #Function to map the covariate
  DFSim <- data.frame(x=grid[,1], y=grid[,2])
  Const.cov<-nndistance
  
  ggplotData <- ggplot(data=sampData, mapping=aes(x=x,y=y)) + 
    geom_tile(data=DFSim, mapping=aes(x=x,y=y,fill=u_ver)) +
    scale_fill_viridis_c(option = "B",direction = -1) +
    geom_point(aes(size = Const.cov),shape = 21,colour = "black", fill = "white") +
    theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
    labs(fill="S(x)")
  ggplotData
}



plot.simulated.field<-function(grid,u_ver){
  #This function plots the simulated S(x)
  
  ggplot() + geom_tile(data=data.frame(x=grid[,1], y=grid[,2],u=u_ver), 
                       mapping=aes(x=x,y=y,fill=u)) +
    scale_fill_viridis_c(option = "B",direction = -1) + theme_bw() +
    labs(fill="S(x)")
}


fit.simple.model<-function(range0,pr=0.01,sigma0,ps=0.01,mesh,sampData){
  #This function fits a single geostatistical model to data 
  
  x<-sampData[,1]
  y<-sampData[,2]
  marks<-sampData[,3]
  
  spde<-inla.spde2.pcmatern(mesh,alpha=2,prior.range=c(range0,pr),prior.sigma=c(sigma0,ps))
  
  formula1 <- y ~ -1 + Intercept + f(spatial.field, model=spde)
  
  s.index1<-inla.spde.make.index(name="spatial.field",n.spde=spde$n.spde)
  
  # Projector matrix from rough mesh for fast computation
  proj_dados1 <- inla.spde.make.A(mesh = mesh, loc =cbind(x,y))
  
  
  ## Stack
  stack.est1 <- inla.stack(data=list(y=marks),
                         A=list(proj_dados1,1),
                         effects=list(i=1:spde$n.spde,
                                      Intercept=rep(1,length(x))),
                         tag='est')
  
  formula.1 <- y ~ -1 + Intercept + f(i, model = spde)
  res1<-inla(formula.1, data=inla.stack.data(stack.est1), 
           control.predictor=list(A=inla.stack.A(stack.est1)),
           control.compute=list(dic=T,waic=T,config=T, cpo=T))
  
  return(res1)
}


  
fit.construct.model<-function(range0,pr=0.01,sigma0,ps=0.01,mesh,sampData,nndistance){
  #This function fits a single model to data including the constructed covariate
  #range0<-0.05  ##VERRR
  #pr<-0.01
  #sigma0<-1
  #ps<-0.01
  
  x<-sampData[,1]
  y<-sampData[,2]
  marks<-sampData[,3]
  
  spde<-inla.spde2.pcmatern(mesh,alpha=2,prior.range=c(range0,pr),prior.sigma=c(sigma0,ps))
  
  #Projector matrix
  proj_dados1 <- inla.spde.make.A(mesh = mesh, loc =cbind(x,y))
  
  
  #Create stack data for inla
  s.index<-inla.spde.make.index(name="spatial.field",n.spde=spde$n.spde)
  
  stack_smooth_dados2 <- inla.stack(data=data.frame(y=marks),
                                    A=list(2,proj_dados1),
                                    effects=list(data.frame(Intercept=rep(1,length(x)),
                                                            cov1=nndistance),
                                                 spatial.field=s.index),
                                    tag='obs')
  
  
  stack_dados2 <- inla.stack(stack_smooth_dados2)
  formula_smooth_dados2 <- y ~ -1 + Intercept + cov1 + f(spatial.field, model=spde)
  
  result_smooth_dados2 <- inla(formula_smooth_dados2,
                               data=inla.stack.data(stack_dados2, spde = spde),
                               family="gaussian",
                               control.predictor = list(A = inla.stack.A(stack_dados2),
                                                        compute = TRUE),
                               num.threads = 2,
                               control.compute=list(dic=T,waic=T,config=T, cpo=T))
                               #,verbose = T)
  
  aux1 <-result_smooth_dados2$summary.fitted.values[inla.stack.index(stack_dados2,"obs")$data, "mean"]
  
  aux2 <-result_smooth_dados2$summary.fitted.values[inla.stack.index(stack_dados2,"obs")$data, "sd"]
  
  resid1 <- (marks - aux1) /aux2
  
  return(list(result_smooth_dados2=result_smooth_dados2,resid1=resid1))
  
  
}


book.mesh.dual <- function(mesh) {
  #From SPDE Krainski Et al. book 
  if (mesh$manifold=='R2') {
    ce <- t(sapply(1:nrow(mesh$graph$tv), function(i)
      colMeans(mesh$loc[mesh$graph$tv[i, ], 1:2])))
    library(parallel)
    pls <- mclapply(1:mesh$n, function(i) {
      p <- unique(Reduce('rbind', lapply(1:3, function(k) {
        j <- which(mesh$graph$tv[,k]==i)
        if (length(j)>0) 
          return(rbind(ce[j, , drop=FALSE],
                       cbind(mesh$loc[mesh$graph$tv[j, k], 1] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 1], 
                             mesh$loc[mesh$graph$tv[j, k], 2] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 2])/2))
        else return(ce[j, , drop=FALSE])
      })))
      j1 <- which(mesh$segm$bnd$idx[,1]==i)
      j2 <- which(mesh$segm$bnd$idx[,2]==i)
      if ((length(j1)>0) | (length(j2)>0)) {
        p <- unique(rbind(mesh$loc[i, 1:2], p,
                          mesh$loc[mesh$segm$bnd$idx[j1, 1], 1:2]/2 +
                            mesh$loc[mesh$segm$bnd$idx[j1, 2], 1:2]/2, 
                          mesh$loc[mesh$segm$bnd$idx[j2, 1], 1:2]/2 +
                            mesh$loc[mesh$segm$bnd$idx[j2, 2], 1:2]/2))
        yy <- p[,2]-mean(p[,2])/2-mesh$loc[i, 2]/2
        xx <- p[,1]-mean(p[,1])/2-mesh$loc[i, 1]/2
      }
      else {
        yy <- p[,2]-mesh$loc[i, 2]
        xx <- p[,1]-mesh$loc[i, 1]
      }
      Polygon(p[order(atan2(yy,xx)), ])
    })
    return(SpatialPolygons(lapply(1:mesh$n, function(i)
      Polygons(list(pls[[i]]), i))))
  }
  else stop("It only works for R2!")
}

fit.preferential.model<-function(range0,pr=0.01,sigma0,ps=0.01,mesh,sampData){
  #This function fits a preferential model to data 
  #range0<-0.05
  #pr<-0.01
  #sigma0<-1
  #ps<-0.01
  
  
  x<-sampData[,1]
  y<-sampData[,2]
  marks<-sampData[,3]  
  
  spde<-inla.spde2.pcmatern(mesh,alpha=2,prior.range=c(range0,pr),prior.sigma=c(sigma0,ps))
  
  #Projector matrix
  proj_dados1 <- inla.spde.make.A(mesh = mesh, loc =cbind(x,y))


  n<-length(x)
  nv<-mesh$n
  y.pp <- rep(0:1, c(nv, n)) #0 for mesh node and 1 for observation
  
  #Dual mesh
  dmesh <- book.mesh.dual(mesh)

  loc.d<- cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0))

  #The domain polygon can be converted into a SpatialPolygons class as follows:
  domain.polys <- Polygons(list(Polygon(loc.d)), '0')
  domainSP <- SpatialPolygons(list(domain.polys))

  #library(rgeos)
  w <- sapply(1:length(dmesh), function(i) {
    if (gIntersects(dmesh[i, ], domainSP))
      return(gArea(gIntersection(dmesh[i, ], domainSP)))
    else return(0)
  })

  e.pp <- c(w, rep(0, n))

  imat <- Diagonal(nv, rep(1, nv))
  lmat <- inla.spde.make.A(mesh, cbind(x,y))

  A.pp <- rbind(imat, lmat)

  stk2.y <- inla.stack(
    data = list(y = cbind(y=marks, NA), e = rep(0, length(x))),
    A = list(proj_dados1, 1),
    effects = list(i = 1:nv, b0.y = rep(1, length(x))),
    tag = 'resp2')

  stk2.pp <- inla.stack(data = list(y = cbind(NA, y.pp), e = e.pp),
                      A = list(A.pp, 1),
                      effects = list(j = 1:nv, b0.pp = rep(1, nv + n)),
                      tag = 'pp2')

  j.stk <- inla.stack(stk2.y, stk2.pp)

  # Gaussian prior
  gaus.prior <- list(prior = 'gaussian', param = c(0, 2))

  # Model formula
  jform <- y ~ 0 + b0.pp + b0.y + f(i, model = spde) +
    f(j, copy = 'i', fixed = FALSE,
      hyper = list(theta = gaus.prior))
  # Fit model
  j.res <- inla(jform, family = c('gaussian', 'poisson'),
              data = inla.stack.data(j.stk),
              E = inla.stack.data(j.stk)$e,
              control.compute=list(dic=T,waic=T,config=T, cpo=T),
              control.predictor = list(A = inla.stack.A(j.stk)))
  
  return(j.res)
}


paper.simul.square1<-function(total_rep=100,a=0.3,b=0.15,sigma0=sqrt(1.5),range0=0.15,mu.field=4,
                             nugget.prec=1/0.1,n=100,beta=2,perc=5){
  #Function to simulate spatial data preferentially (beta is the preferential parameter) in a 
  #square domain [0,1]x[0,1], test whether that preferentiability is spotted (MLC test), to
  #construct a covariate to mitigate the preferential sampling effect and to test whether it
  #accomplished its purposes when included in modelling, by testing the residuals of the model 
  #that includes that covariate for preferentiability 
  
  #Inputs
  #Parameters used to construct the mesh to simulate the field
  #total_rep - number of replicas to do
  #a - inner max.edge
  #b - outer max.edge
  #sigma0 - marginal sd of the field
  #range0 - scale
  #mu.field - mean of the random field S(x)
  #beta - preferential parameter 
  #nugget.prec - nugget precision (of Y)
  #n - number of observations to sample
  
  
  
  # mesh to simulate spatially structured random effect
  # the study region is set as a square domain: (0,1)x(0,1) \subset R^2
  domain_boundary <- matrix(c(0,0,1,0,1,1,0,1,0,0), byrow=TRUE, ncol=2) 
  mesh <- inla.mesh.2d(loc.domain=domain_boundary, max.edge=c(a, b))
  #Plot the mesh - use function gg of library(inlabru)
  ggplot() + gg(mesh) + theme(axis.title.x = element_blank(),axis.title.y = element_blank())
  #VER!!! gg()
  
  rep<-0
  while(rep < total_rep){
    rep=rep+1
    
    #Simulation of the spatial random effect of mean mu.field, in a given grid in the unit square domain
    sim.ran.eff<-simul.rand.effect(range0,sigma0,mesh=mesh,mu.field)
    u_ver<-sim.ran.eff$u_ver
    grid.REff<-sim.ran.eff$grid  #Grid of the simulated random effect
    
    
    #Plot simulated random field
    p1<-plot.simulated.field(grid.REff,u_ver)
    p1
    
    #Sample the data according to PP with intensity exp(beta*S(x))
    sampData<-sample.prefer(u_ver,nugget.prec,grid.REff,n,beta)
    
    #Plot the sample data over the surface S(x)
    p1+geom_point(data=sampData,shape = 21,colour = "black", fill = "white", size=1,mapping=aes(x=x,y=y)) +
      theme(axis.title.x = element_blank(),axis.title.y = element_blank())
    
    #MLC test for data preferentiability
    RES<-MLC(domain_boundary,sampData)
    
    #Write results out
    write.table(matrix(RES,ncol=3), 
                file=paste("RES_nocovariate_beta",beta,"_",n,"obs_",perc,"p.txt",sep=""), 
                append=T,row.names=F, col.names=F)
    
    
    
    #Obtain the construct covariate Distance to the k nearest neighbours (5% nearest neighbours)
    nndistance<-construct.covariate(sampData,perc)
    
    
    #Plot the covariate
    plot.covariate(sampData,nndistance,grid.REff,u_ver)
    
    #Fit model to data including the constructed covariate, no adjustment for preferential sampling
    res.aux<-fit.construct.model(range0,pr=0.01,sigma0,ps=0.01,mesh,sampData,nndistance)
    result_smooth_dados2<-res.aux$result_smooth_dados2
    
    summary(result_smooth_dados2)
    
    resid1<-data.frame(x=sampData$x,y=sampData$y,resid=res.aux$resid)  #A coluna 3 do sampData não deveria chamar-se ysim; Impacta o ficheiro do teste
    
    #MLC test for resuduals preferentiability
    RES2<-MLC(domain_boundary,resid1)
    
    #Write results out
    write.table(matrix(RES2,ncol=3), 
                file=paste("RES_covariate_beta",beta,"_",n,"obs_",perc,"p.txt",sep=""), 
                append=T,row.names=F, col.names=F)
    
  } #end while
  
  
}



paper.simul.square2<-function(total_rep=100,a=0.3,b=0.15,sigma0=sqrt(1.5),range0=0.15,mu.field=4,
                              nugget.prec=1/0.1,n=100,beta=2,perc=5){
  #Function to simulate spatial data preferentially (beta is the preferential parameter) in a 
  #square domain [0,1]x[0,1], test whether that preferentiability is spotted (MLC test), 
  #to fit a simple geostatistical model to simulated data
  #to
  #construct a covariate to mitigate the preferential sampling effect and to test whether it
  #accomplished its purposes when included in modelling, by testing the residuals of the model 
  #that includes that covariate for preferentiability 
  
  #Inputs
  #Parameters used to construct the mesh to simulate the field
  #total_rep - number of replicas to do
  #a - inner max.edge
  #b - outer max.edge
  #sigma0 - marginal sd of the field
  #range0 - scale
  #mu.field - mean of the random field S(x)
  #beta - preferential parameter 
  #nugget.prec - nugget precision (of Y)
  #n - number of observations to sample
  
  
  
  # mesh to simulate spatially structured random effect
  # the study region is set as a square domain: (0,1)x(0,1) \subset R^2
  domain_boundary <- matrix(c(0,0,1,0,1,1,0,1,0,0), byrow=TRUE, ncol=2) 
  mesh <- inla.mesh.2d(loc.domain=domain_boundary, max.edge=c(a, b))
  #Plot the mesh - use function gg of library(inlabru)
  ggplot() + gg(mesh) +
    theme(axis.title.x = element_blank(),axis.title.y = element_blank())    
  ###TENTAR TIRAR DEPENDENCIA DO GG!!!
  
  
  rep<-0
  while(rep < total_rep){
    rep=rep+1
    
    #Simulation of the spatial random effect of mean mu.field, in a given grid in the unit square domain
    sim.ran.eff<-simul.rand.effect(range0,sigma0,mesh=mesh,mu.field)
    u_ver<-sim.ran.eff$u_ver
    grid.REff<-sim.ran.eff$grid  #Grid of the simulated random effect
    
    
    ##Plot simulated random field
    #p1<-plot.simulated.field(grid.REff,u_ver)
    #p1
    
    #Sample the data according to PP with intensity exp(beta*S(x))
    sampData<-sample.prefer(u_ver,nugget.prec,grid.REff,n,beta)
    
    ##Plot the sample data over the surface S(x)
    #p1+geom_point(data=sampData,shape = 21,colour = "black", fill = "white", size=1,mapping=aes(x=x,y=y)) +
    # theme(axis.title.x = element_blank(),axis.title.y = element_blank())
    
    #MLC test for data preferentiability
    RES<-MLC(domain_boundary,sampData)
    
    #Write results out
    write.table(matrix(RES,ncol=3), 
                file=paste("RES_nocovariate_beta",beta,"_",n,"obs_",perc,"p.txt",sep=""), 
                append=T,row.names=F, col.names=F)
    
    #To fit a simple geostatistical model
    res1<-fit.simple.model(range0,pr=0.01,sigma0,ps=0.01,mesh,sampData)
    
    
    resid1<-data.frame(x=sampData$x,y=sampData$y,resid=res1$resid)  #A coluna 3 do sampData não deveria chamar-se ysim; Impacta o ficheiro do teste
    
    #MLC test for residuals preferentiability
    RES1<-MLC(domain_boundary,resid1)
    
    #Write results out
    write.table(matrix(RES1,ncol=3), 
                file=paste("RES_simple_model_beta",beta,"_",n,"obs_",perc,"p.txt",sep=""), 
                append=T,row.names=F, col.names=F)
    
    
    #Obtain the construct covariate Distance to the k nearest neighbours (5% nearest neighbours)
    nndistance<-construct.covariate(sampData,perc)
    
    
    ##Plot the covariate
    #plot.covariate(sampData,nndistance,grid.REff,u_ver)
    
    #Fit model to data including the constructed covariate, no adjustment for preferential sampling
    res.aux<-fit.construct.model(range0,pr=0.01,sigma0,ps=0.01,mesh,sampData,nndistance)
    result_smooth_dados2<-res.aux$result_smooth_dados2
    
    summary(result_smooth_dados2)
    
    #Quality Criteria
    cpo2<--mean(log(result_smooth_dados2$cpo$cpo), na.rm=T)
    
    write.table(matrix(c(cpo2),ncol=1), 
                file=paste("criteria_",n,"obs_beta_",beta,"covariate_model.txt",sep=""), append=T,row.names=F, col.names=F)
    
    
    resid2<-data.frame(x=sampData$x,y=sampData$y,resid=res.aux$resid)  #A coluna 3 do sampData não deveria chamar-se ysim; Impacta o ficheiro do teste
    
    #MLC test for residuals preferentiability
    RES2<-MLC(domain_boundary,resid2)
    
    #Write results out
    write.table(matrix(RES2,ncol=3), 
                file=paste("RES_covariate_beta",beta,"_",n,"obs_",perc,"p.txt",sep=""), 
                append=T,row.names=F, col.names=F)
    
    
    #Fit a preferential model
    j.res<-fit.preferential.model(range0,pr=0.01,sigma0,ps=0.01,mesh,sampData)
    summary(j.res)
    
    #Quality Criteria
    cpo3<--mean(log(j.res$cpo$cpo), na.rm=T)
    
    write.table(matrix(c(cpo3),ncol=1), 
                file=paste("criteria_",n,"obs_beta_",beta,"preferential_model.txt",sep=""), append=T,row.names=F, col.names=F)
    
    
  } #end while
  
  
}

