# population structure: (S, E, E_t, I, I_t, I_s, H, R, D)
#   S = susceptible
#   E = latent
#   E_t = latent but test-posI_tive
#   I = infectious
#   I_t = infectious and test-posI_tive
#   I_s = symptomatic
#   H = hospI_talized
#   R = recovered
#   D = dead
#####################
# Assumptions:
# 1. latent and infectious individuals can be identified by testing
#   1a. test-positive infectious individuals will reduce contacts 90% (modifiable)
# 2. latent individuals will become infectious, then symptomatic before seeking hospital care
#   2a. only a portion (age-specific) of symptomatic individuals will require hospitalization
# 3. recovery I_s based on age and disease status
#   3a. a portion of latent individuals will recover before being symptomatic
#   3b. symptomatic individuals who are not hospitalized will recover
#   3c. hospitalized individuals will either recover or die

params=list(
  # most from http://gabgoh.github.io/COVID/index.html
  l=5,   #length of latency
  p_s=0.96, #probability of symptoms
  i=3,   #length of infectiousness before symptomatic
  r=11,   #time to recovery if mild/asymptomatic
  p_h=0.14, #probability of hospitalization given symptoms, from 10.1542/peds.2020-0702
  s=5,   #length of symptoms before hospitalization
  h=10,  #length of hospitalization
  p_d=0.04 #probability of death given hospitalization
)
paramlist=list( #set up separate parameters for each age group
  params,params,params,params,params,params,params,params,params
)
paramlist[[1]]$p_s=0.1 #from 10.1542/peds.2020-0702
#other parameters can be made age-specific as needed

lambda_params=list(
  # from http://gabgoh.github.io/COVID/index.html
  p_i=.45,  #probability of infection given contact
  p_hcw=.45, #probability of infection of HCW given contact (could change after PPE decrease)
  c_hcw=45, #number of contacts by HCW per hospitalized patient per day
  # modifiable 
  q=0.9, #proportional decrease in c due to quarantine after positive test
  d=0, #proportional decrease in c due to social distancing
  w=0 #proportional decrease in p_i due to hygiene
)

# from https://journals.plos.org/plosmedicine/article/file?id=10.1371/journal.pmed.0050074&type=printable
waifw=read.csv("contact_matrix.csv")

cfr_byage=data.frame(
  #from https://www.cdc.gov/mmwr/volumes/69/wr/mm6912e2.htm
  #assume HCW similar to 22-54
  agegroup=colnames(waifw),
  p_h=c(2,.175,.175,.2475,.253,.36,.446,.508,.21),
  cfr=c(0,.0015,.0015,.0065,.02,.038,.074,.19,.004)
)
cfr_byage$p_d=cfr_byage$cfr/cfr_byage$p_h #convert p(d) to p(d|h)

#pops is a matrix with the current population in each age category
#(S, E, E_t, I, I_t, I_s, H, R, D)

Lambda=function(lambda_params,waifw,pops){
  lambda=rep(NA,nrow(pops))
  beta=waifw*lambda_params$p_i
  I_mat=apply(cbind(pops[1:8,4]*(1-lambda_params$d), #untested infectious
                      pops[1:8,5:6]*(1-lambda_params$q)), #quarantined
              1,sum)/apply(pops[1:8,1:8],1,sum) #add up pressure from infectious groups, normalizes by group size
  for(i in 1:ncol(pops)){
    lambda[i]=sum(beta[1:8,i]*I_mat) #sum up infectious pressure from each age group
  }
  lambda[9]=lambda_params$p_hcw*sum(pops[,7])/sum(pops[9,])*lambda_params$c_hcw #average all hospitalized patients across HCWs
  return(lambda)
}


DZstep=function(pop,params,lambda){
  S=pop[1]; E=pop[2]; E_t=pop[3]; I=pop[4]; I_t=pop[5]; I_s=pop[6]; H=pop[7]; R=pop[8]; D=pop[9]
  # hospitalized
    newdeaths=min(H,rpois(1,params$p_d*(1/params$h)*H)) #how many will die
      D=D+newdeaths;H=H-newdeaths
    recoverH=min(H,rpois(1,(1-params$p_d)*(1/params$h)*H)) #hospitalized recover
      R=R+recoverH;H=H-recoverH
  # symptomatic
    hospitalize=min(I_s,rpois(1,params$p_h*(1/params$s)*I_s)) #symptomatic become hospitalized
      H=H+hospitalize;I_s=I_s-hospitalize
    recoverI_s=min(I_s,rpois(1,(1-params$p_h)*(1/params$r)*I_s)) #symptomatic recover at home
      R=R+recoverI_s;I_s=I_s-recoverI_s
  # infectious
    symptomsI=min(I,rpois(1,I*(1/params$i)*params$p_s)) #infectious become symptomatic
    symptomsI_t=min(I_t,rpois(1,I_t*(1/params$i)*params$p_s)) #tested infectious become symptomatic
      I_s=I_s+symptomsI+symptomsI_t;I=I-symptomsI;I_t=I_t-symptomsI_t
    recoverI=min(I,rpois(1,I*(1/params$i)*(1-params$p_s))) #infectious recover
    recoverI_t=min(I_t,rpois(1,I_t*(1/params$i)*(1-params$p_s))) #tested infectious recover
      R=R+recoverI+recoverI_t;I=I-recoverI;I_t=I_t-recoverI_t
  # latent
    infectious=min(E,rpois(1,E/params$l)) #latent become infectious
      I=I+infectious;E=E-infectious
    infectious_t=min(E_t,rpois(1,E_t/params$l)) #tested latent become infectious
      I_t=I_t+infectious_t;E_t=E_t-infectious_t
  # susceptible
    infection=min(S,rpois(1,S*lambda)) #susceptible become infected
      S=S-infection;E=E+infection
  #recover population
    newpop=c(S,E,E_t,I,I_t,I_s,H,R,D)
    return(list(pop=newpop,deaths=newdeaths,infections=infection))
}