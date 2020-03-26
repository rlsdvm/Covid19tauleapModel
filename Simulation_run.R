#####Running a simulation
source("SEIR_functions.R")

###############No control, Chicago
N=2705994 #total population of Chicago city per census
agedist=c(.22,.09,.32,.12,.11,.08,.04,.02, #proportion in each age category, per census
        175100) #number of physicians, RN, NA, MA, and LPN/LVN, per https://matter.health/uploads/hc3_whitepaper_digital.pdf
agenums=c(round((N-agedist[9])*agedist[1:8],0),agedist[9]) #distribute numbers across age groups

timestop=30*6 #run for 6 months

####one iteration
RunIter=function(N,agenums,timestop,waifw,params,lambda_params){
poparray=array(NA,dim=c(9,9,timestop),dimnames = list(colnames(waifw),c("S","E","E_t","I","I_t","I_s","H","R","D"),paste("t",1:timestop,sep="")))
    #sets up an array for the population at each timestep in each age and disease category
poparray[,,1]=0
  poparray[,1,1]=agenums #set up the starting population as fully susceptible
  startdz=100 #how many diseased to start with
  startdist=rmultinom(1,startdz,agedist[4:8]) #distribute the diseased across the older age categories
  poparray[4:8,1,1]=poparray[4:8,1,1]-startdist #take diseased out of S
  poparray[4:8,4,1]=startdist #put diseased in I
 #initialize saving
  deaths=infections=lambda=matrix(NA,nrow=timestop,ncol=9)
    colnames(deaths)=colnames(infections)=colnames(lambda)=colnames(waifw)
    seir=matrix(NA,nrow=timestop,ncol=9);colnames(seir)=c("S","E","E_t","I","I_t","I_s","H","R","D")
 #run the simulation
  for(t in 2:timestop){ #daily
    lambda[t,]=Lambda(lambda_params,waifw,poparray[,,t-1])
    #step each agegroup through infections
    for(g in 1:9){
      tstep=DZstep(poparray[g,,t-1],paramlist[[g]],lambda[t,g])
      deaths[t,g]=tstep$deaths
      infections[t,g]=tstep$infections
      poparray[g,,t]=tstep$pop
    }
  seir[t,]=apply(poparray[,,t],2,sum)
    }
    return(list(poparray=poparray,seir=seir,deaths=deaths,infections=infections))
}

#######multiple iterations
n_iter=1000 #number of iterations
iter_out=list(NULL) #initialize
seir_array=array(NA,dim=c(timestop,9,n_iter),
                 dimnames = list(paste("t",1:timestop,sep=""),c("S","E","E_t","I","I_t","I_s","H","R","D"),1:n_iter))
death_array=infection_array=array(NA,dim=c(timestop,9,n_iter),
                dimnames = list(paste("t",1:timestop,sep=""),colnames(waifw),1:n_iter))

for(i in 1:n_iter){
  iter_out[[i]]=RunIter(N,agenums,timestop,waifw,params,lambda_params)
  seir_array[,,i]=iter_out[[i]]$seir
  death_array[,,i]=iter_out[[i]]$deaths
  infection_array[,,i]=iter_out[[i]]$infections
}

matplot(infection_array[,7,],type = "l",col = "gray") #plot hospitalized
matplot(seir_array[,9,],type = "l",col = "gray") #plot deaths

    