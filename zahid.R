################################################################################
### Simplified global carbon cycle  (SGCC) model
##################################################################################
### Code superviser: Andreas Will, will@b-tu.de
### Code modifier: Md Zahid Hasan,3814489
show(paste("version env.R:",version,date))
#################################################
### Model Usage: ################################
### see "GCC_ESM_AW_config_Vx"
#################################################
### Model Documentation: SGCCM_users-guide_Vx.pdf
### Model History: SGCCM_README_CHANGES.txt
##################################################################################
### MODEL BEGIN
##################################################################################
show(paste("run.R started"))
#################################################################################
### emission scenario loop ######################################################
#################################################################################
debug=0
epsilon=0.000000000000001

for (k in nr1:nr2){
  #while (deltacomp > epsilon ){
  ### scanario loop for optimisation of initial conditions ######################
  ### if ( k == nr1){
  ###############################################################################
  ### INITIALIZATION OF THE FIELDS
  ###############################################################################
  kp=k+1
  diagtab[5,1]=-9999
  diagtab[6,1]=-9999
    iterres=0
    m=0
    par0[1]=na_i
    par0[2]=nb_i
    par0[3]=nu_i
    par0[4]=nd_i
    par0[5]=kab_i
    par0[6]=kba_i
    par0[7]=kau_i
    par0[8]=kua_i
    par0[9]=kud_i
    par0[10]=kdu_i
    deltagcc<-deltagcc_i
    normtab<-normtab_i
    PARAMFILE<-paste("paramtab_",colnames(ghgemis)[k],sep="")
    if (warmstart == 1){
      paramtab<-read.table(paste(RESULTDIR,"/",PARAMFILE,sep=""))
      for (l in np1:np2){
        nl=l-np1+1
        par0[l]=paramtab[length(paramtab[,1]),nl]
      }
    }
    na[1,]=par0[1]
    nb[1,]=par0[2]
    nu[1,]=par0[3]
    nd[1,]=par0[4]
    kab=par0[5]
    kba=par0[6]
    kau=par0[7]
    kua=par0[8]
    kud=par0[9]
    kdu=par0[10]
    ###
    for(i in 1:ntiter){
      na[i+1,k]=na[i,k]+dt*(-(kab+kau)*na[i,k]+kba*nb[i,k]+kua*nu[i,k]+ghgemis[i,k]) 
      nb[i+1,k]=nb[i,k]+dt*((kab*na[i,k])-kba*nb[i,k]) 
      nu[i+1,k]=nu[i,k]+dt*((kau*na[i,k])-(kud+kua)*nu[i,k]+kdu*nd[i,k]) 
      nd[i+1,k]=nd[i,k]+dt*(kud*nu[i,k]-kdu*nd[i,k])
      #if (i == 1 || i == ndt)  show(paste("ndt=",i," l=",l," k=",k," na=",na[i,k]," nb=",nb[i,k]," nu=",nu[i,k]," nd=",nd[i,k]))
      # norm calculation
      diff_na[i] = na[i+1,k]-ghgconcgt[i+1,k]
      error_na[i] = (abs(diff_na[i]))^expon
      if (Lnorm == 3){ error_na[i] = (log(na[i+1,k]/ghgconcgt[i+1,k]))^expon}
      if (Lnorm == 4){ error_na[i] = (abs(diff_a[i]/ghgconcgt[i+1,k]))^expon}
    }
    for (i in 1:ns){
      for (j in 1:np){
        spreadgcc[i,j]<-deltagcc[j]*spread0[i]
      }
    }
    ### iteration loop (iter) #############################################
    for (iter in n1iter:n2iter){
      iterres=iterres+1
      if (iterres == niterincdoc+1){ iterres=1}
      for (l in np1:np2){
        minlist[l]=min(normtab[l])
      }
      nmin<-which.min(minlist[np1:np2])+np1-1
      for (l in np1:np2){
        if (nmin == l){ 
          par0[l] = par0[l]*(1 +spreadgcc[which.min(normtab[,l]),l])
          if (par0[l] < parmin[l]) ( par0[l]=parmin[l])   ### keep minimum parameter value
          if (par0[l] > parmax[l]) ( par0[l]=parmax[l])   ### keep maximum parameter value
          if ( ( which.min(normtab[,l]) == nsmean ||  par0[l] == parmin[l] || par0[l] == parmax[l] ) && ( deltagcc[l] > epsilon ) ) {
            deltagcc[l]=deltagcc[l]/qdelta 
            deltacomp=deltagcc[l]
            if (debug > 0 ) {show(paste("decdelta: iter, l, deltagcc:",iter,l,deltagcc[l]))}
          }
	      }
        if (debug > 0) {
          show(paste("decdelta: iter, l, which.min(normtab[,l]), parmin,par0,parmax:",iter,l, which.min(normtab[,l]),parmin[l],par0[l],parmax[l])) 
	        if ( which.min(normtab[,l]) == 1 ) {
		        show(paste("condition 1"))
	        }
 	        if ( abs(normtab[1,l] - normtab[nsmean,l]) > epsilon ){
		        show(paste("condition 2"))
	        }
   	      if ( which.min(normtab[,l]) == ns ) {
		        show(paste("condition 3"))
	        }
 	        if ( abs(normtab[ns,l] - normtab[nsmean,l]) > epsilon ) {
		        show(paste("condition 4"))
	        }
        }
        #if ( ( which.min(normtab[,l]) == 1 && ( abs(normtab[1,l]- normtab[nsmean,l]) > epsilon ) ) || ( which.min(normtab[,l]) == ns && (abs(normtab[ns,l]- normtab[nsmean,l]) > epsilon ) ) ) {
        # deltagcc[l]=deltagcc[l]*qdelta
        # if (debug > 0 ) {show(paste("incdelta: iter, l, deltagcc:",iter,l,deltagcc[l]))}
        #}
        if (debug > 0 ) {
      		show(paste("incdelta: iter, l, which.min(normtab[,l]), normtab(1,nsmean,ns):",iter,l,which.min(normtab[,l]),normtab[1,l],normtab[nsmean,l],normtab[ns,l])) 
        }
      }
      deltafmax=9.999
      for (l in np1:np2) {
        if (deltagcc[l]/deltafmax > deltacomp) (deltagcc[l]=deltagcc[l]/qdelta)
        for (j in 1:ns){
          spreadgcc[j,l]<-deltagcc[l]*spread0[j]
        }
      }
      ### parameter loop : l ##################################################
      for (l in np1:np2){
        ### spread loop #######################################################
        na[1,k]=par0[1]
        nb[1,k]=par0[2]
        nu[1,k]=par0[3]
        nd[1,k]=par0[4]
        kab=par0[5]
        kba=par0[6]
        kau=par0[7]
        kua=par0[8]
        kud=par0[9]
        kdu=par0[10]
        for (j in 1:ns){
          if (l == 1) { na[1,k]=par0[1]*(1+spreadgcc[j,1]) }
          if (l == 2) { nb[1,k]=par0[2]*(1+spreadgcc[j,2]) }
          if (l == 3) { nu[1,k]=par0[3]*(1+spreadgcc[j,3]) }
          if (l == 4) { nd[1,k]=par0[4]*(1+spreadgcc[j,4]) }
          if (l == 5) { kab=par0[5]*(1+spreadgcc[j,5]) }
          if (l == 6) { kba=par0[6]*(1+spreadgcc[j,6]) }
          if (l == 7) { kau=par0[7]*(1+spreadgcc[j,7]) }
          if (l == 8) { kua=par0[8]*(1+spreadgcc[j,8]) }
          if (l == 9) { kud=par0[9]*(1+spreadgcc[j,9]) }
          if (l == 10){ kdu=par0[10]*(1+spreadgcc[j,10]) }
          ### time loop ###################################################################
          for(i in 1:ntiter){
            na[i+1,k]=na[i,k]+dt*(-(kab+kau)*na[i,k]+kba*nb[i,k]+kua*nu[i,k]+ghgemis[i,k]) 
	          nb[i+1,k]=nb[i,k]+dt*((kab*na[i,k])-kba*nb[i,k]) 
	          nu[i+1,k]=nu[i,k]+dt*((kau*na[i,k])-(kud+kua)*nu[i,k]+kdu*nd[i,k]) 
	          nd[i+1,k]=nd[i,k]+dt*(kud*nu[i,k]-kdu*nd[i,k])
	   # if (i == 1 || i == ndt)  show(paste("ndt=",i," l=",l," k=",k," na=",na[i,k]," nb=",nb[i,k]," nu=",nu[i,k]," nd=",nd[i,k]))
	   # norm calculation
  	        diff_na[i] = na[i+1,k]-ghgconcgt[i+1,k]
	          error_na[i] = (abs(diff_na[i]))^expon
	          if (Lnorm == 3){ error_na[i] = (log(na[i+1,k]/ghgconcgt[i+1,k]))^expon}
	          if (Lnorm == 4){ error_na[i] = (abs(diff_a[i]/ghgconcgt[i+1,k]))^expon}
          }
          ### time loop end ###############################################################
          normtab[j,l]= (sum(error_na[n1:n2])/(n2-n1+1))^(1/expon) 
          (sum(error_na[n1:n2])/(n2-n1+1))^(1/expon)
        }
        ### parameter loop ( l ) end ###################################################
      }
      ### iteration loop (iter) end #############################################
      if (iterres == 1) { 
	      m=m+1
      	paramtab[m,]<-c(par0[np1:np2],normtab[nsmean,np1:np2],deltagcc[np1:np2])
      	show(paste(format(c(k,iter,par0[np1:np2],normtab[nsmean,np1:np2],deltagcc[np1:np2]),digits=4)))
      	if ( warmstart == 0 || ( warmstart == 1 && n2iter > n1iter )){
      	  if (debug > 0) show(paste("FILE=",PARAMFILE))
      	  write.table(format(paramtab,digits=4), file=paste(RESULTDIR,"/",PARAMFILE,sep=""))
      	}
      }
      ### scanario loop for optimisation of initial conditions end ######################
    }
  ###}
  ### Generate the time series for optimum parameter values:
  na[1,k]=par0[1]
  nb[1,k]=par0[2]
  nu[1,k]=par0[3]
  nd[1,k]=par0[4]
  kab=par0[5]
  kba=par0[6]
  kau=par0[7]
  kua=par0[8]
  kud=par0[9]
  kdu=par0[10]
  resulttab[1,1]=na[1,1]
  resulttab[1,2]=nb[1,1]
  resulttab[1,3]=nu[1,1]
  resulttab[1,4]=nd[1,1]
  ### time loop ###################################################################
  carbonm[1,k]=par0[1]+par0[2]+par0[3]+par0[4]
  carbonbudget[1,k]=carbonm[1,k]-ghgemis[1,k]-carbonm[1,k]
  ghgemisac[1,k]=ghgemis[1,k]
  for(i in 1:ndt){
    na[i+1,k]=na[i,k]+dt*(-(kab+kau)*na[i,k]+kba*nb[i,k]+kua*nu[i,k]+ghgemis[i,k]) 
    nb[i+1,k]=nb[i,k]+dt*((kab*na[i,k])-kba*nb[i,k]) 
    nu[i+1,k]=nu[i,k]+dt*((kau*na[i,k])-(kud+kua)*nu[i,k]+kdu*nd[i,k]) 
    nd[i+1,k]=nd[i,k]+dt*(kud*nu[i,k]-kdu*nd[i,k])
    resulttab[i+1,1]=na[i+1,k]
    resulttab[i+1,2]=nb[i+1,k]
    resulttab[i+1,3]=nu[i+1,k]
    resulttab[i+1,4]=nd[i+1,k]
    #   show(paste("i, na,nb,nu,nd",i,na[i,k],nb[i,k],nu[i,k],nd[i,k]))
    carbonm[i+1,k]=na[i+1,k]+nb[i+1,k]+nu[i+1,k]+nd[i+1,k]
    ghgemisac[i+1,k]=ghgemisac[i,k]+ghgemis[i+1,k]
    carbonbudget[i+1,k]=carbonm[i+1,k]-ghgemisac[i,k]-carbonm[1,k]
    if (debug > 0) show(paste(na[i+1,k],diagtab[6,kp]))
    if (na[i+1,k] > diagtab[6,kp] ) {
      if (debug > 0) show(paste("diag: i, k, na:",i,k,na[i,k]))
      diagtab[5,kp]=i*dt+tstart
      diagtab[6,kp]=na[i+1,k]
    }
    if (debug > 0) show(paste("diag: k, t, na(t), tmax, na(tmax), tref, na(tref): ",k,tstart+i*dt, na[i,k], diagtab[5,kp],diagtab[6,kp],diagtab[7,kp],diagtab[8,kp]))
  }
  diagtab[1,kp]=na[1,k]
  diagtab[2,kp]=na[ndt+1,k]
  diagtab[3,kp]=na[t1-tstart+1,k]
  diagtab[4,kp]=na[t2-tstart+1,k]
  for (ntr in 1:length(tref)){
    ndiag=length(diagr1)+ntr
    diagtab[ndiag,kp]=na[tref[ntr]-tstart,k]
  }
  format(par0[1], scientific=FALSE, digits= 2)
  DIAGFILE<-paste("diagtab_",
                  format(par0[1], scientific=FALSE, digits= 2),"-",
                  format(par0[2], scientific=FALSE, digits= 2),"-",
                  format(par0[3], scientific=FALSE, digits= 2),"-",
                  format(par0[4], scientific=FALSE, digits= 2),"-",
                  format(par0[5], scientific=FALSE, digits= 2),"-",
                  format(par0[6], scientific=FALSE, digits= 2),"-",
                  format(par0[7], scientific=FALSE, digits= 2),"-",
                  format(par0[8], scientific=FALSE, digits= 2),"-",
                  format(par0[9], scientific=FALSE, digits= 2),"-",
                  format(par0[10], scientific=FALSE, digits= 2),
                  sep="")
  ### time loop end #######################################################
  RESULTFILE<-paste("resulttab_",colE[k],"_",tstart,"-",tend, sep="")
  write.table(format(resulttab,digits=6), file=paste(RESULTDIR,"/",RESULTFILE,sep=""))
}
### emission scenario loop end ############################################
write.table(format(diagtab,digits=4), file=paste(RESULTDIR,"/",DIAGFILE,sep=""))
show(paste( "run.R finished, number of iterations: ",n2iter))
##################################################################
#### run finshed #################################################
##################################################################
