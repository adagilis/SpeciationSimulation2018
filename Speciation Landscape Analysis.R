#Author: Andrius J. Dagilis
#Date: August 2018
#Email questions to adagilis@gmail.com

###Introduction
#These are the scripts used to generate the simulations and figures for "The evolution of hybrid fitness during speciation"
#Some code is optimized, but in general not. Simulations and analyses require a large amount of RAM as a result.
#Not all the results produced here were used in the paper. Likewise, some parts of the script were
#rerun using different parameters, and these runs are not annotated here. 


#In your working directory you should have all of the files from the cellmap.org complete data-set
#Libraries
library(ggplot2)
library(cowplot)
library(wesanderson)
library(reshape2)
library(compiler)
library(foreach)
library(doSNOW)
library(Rmpfr)

#Load Data sets

data.DAmP = read.table("SGA_DAmP.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE,fill=TRUE,quote="")
data.NxN = read.table("SGA_NxN.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE,fill=TRUE,quote="")
data.ExN_NxE = read.table("SGA_ExN_NxE.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE,fill=TRUE,quote="")
data.ExE = read.table("SGA_ExE.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE,fill=TRUE,quote="")

#Read in Ensembl annotated data file. Has plenty of useful info for other analyses (chrom location, etc.)

BioData = read.table("subYeast.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE,fill=TRUE)

#Read in Single knockout fitness file.  

SFit = read.table("strain_ids_and_single_mutant_fitness.csv",header=TRUE,sep=",",stringsAsFactors = FALSE,fill=TRUE,quote="")
#Some strange extra columns that are empty, so get rid of them
SFit = SFit[,1:7]

#For initial analysis, let's only use complete data:
FitData = SFit[complete.cases(SFit),]

#Some universal variables

data.tot = rbind(data.ExE,data.NxN,data.DAmP,data.ExN_NxE)
rm(data.DAmP,data.ExE,data.ExN_NxE,data.NxN)

data.length = length(data.tot[,1])

uniLoc = unique(SFit$Strain.ID)
nLoci = length(uniLoc)

uniStrains = unique(c(data.tot$Query.Strain.ID,data.tot$Array.Strain.ID))
nStrains = length(uniStrains)

strainFits = uniStrains

#The fitnesses in SFit end up not being consistent with those listed in data.tot, so we calculate a mean fitness per strain:
#Defining Strain Fitnesses

idx1 = which(uniStrains %in% data.tot$Query.Strain.ID)
idx2 = which(uniStrains %in% data.tot$Array.Strain.ID)
#Just check for duplicates real quick
idx3 = which(duplicated(c(idx1,idx2)))
#None, moving along, get first occurence of each strain.
idx3 =  !duplicated(data.tot$Query.Strain.ID)
idx4 = !duplicated(data.tot$Array.Strain.ID)
#Match them for fitness:
strainFits[idx1] = data.tot$Query.single.mutant.fitness..SMF.[idx3]
strainFits[idx2] = data.tot$Array.SMF[idx4]
strainFits = sapply(strainFits,as.numeric)

rm(idx1,idx2,idx3,idx4)

Mutations =data.frame(id = uniStrains,s=strainFits-1)
Mutations = Mutations[complete.cases(Mutations),]
nStrains = length(Mutations$id)

#Our DMI metrics, generating complete data set that includes multiple DMI metrics, calculates fixation probability
data = data.tot

#Expected Fitnesses First
ExpFitMult = data$Query.single.mutant.fitness..SMF.*data$Array.SMF
ExpFitAdd = (data$Query.single.mutant.fitness..SMF.+data$Array.SMF-1)
#Calculate Epsilon - there will be some non-sensical values here, but fitness is measured via growth, so there are negative values for fitness. (about 520 of the 20million)
DmiMult = data$Double.mutant.fitness/ExpFitMult-1
DmiAdd = data$Double.mutant.fitness-ExpFitAdd

#Setting extremely negative values to lethal
data$DMIMult[which(data$DMIMult < -1)] = -1
DmiMult[which(DmiMult< -1)] = -1
DmiAdd[which(DmiAdd< -1)] = -1

#Save to data frame
data$DMIMult = DmiMult
data$DMIAdd = DmiAdd
data.tot = data
rm(data,DmiMult,DmiAdd)


#Quality-control: how does epistatic effect relate to p-value reported in Costanzo

ggplot(data,aes(x=P.value,y=DMIMult,z=DmiMult))+stat_summary_2d(fun="length")
#shows that magnitude of interactions is inverse to p-value, so our largest interactions we trust, and minor interactions we're ok with behaving badly.
idx = which(data$P.value !=1)
data = data[idx,]
DmiAdd = DmiAdd[idx]
DmiMult = DmiMult[idx]


#EpsilonTable: A table where every entry is the epistatic interaction between the two strains. Makes look-up faster later
EpsiTable = matrix(nrow=nStrains,ncol=nStrains,0)
EpsiTableAdd = matrix(nrow=nStrains,ncol=nStrains,0)
EpsiTableC = matrix(nrow=nStrains,ncol=nStrains,0)
for(i in 1:(nStrains-1)){
  idx1 = which(data.tot$Query.Strain.ID == Mutations$id[i])
  idx2 = which(data.tot$Array.Strain.ID == Mutations$id[i])
  for(j in (i+1):nStrains){
    idx3 = which(data.tot$Query.Strain.ID[idx2]==Mutations$id[j])
    idx4 = which(data.tot$Array.Strain.ID[idx1]==Mutations$id[j])
    idx = union(idx1[idx4],idx2[idx3])
    if(length(idx)>0){
      EpsiTable[min(i,j),max(i,j)]=mean(data.tot$DMIMult[idx],na.rm=TRUE)
      EpsiTableAdd[min(i,j),max(i,j)] = mean(data.tot$DMIAdd[idx],na.rm=TRUE)
      EpsiTableC[min(i,j),max(i,j)] = mean(data.tot$Genetic.interaction.score..ε.[idx],na.rm=TRUE)
    }
  }
}



#Some important parameters
#Effective population size, will impact drift probability for allele fixation
Ne = 10e06
#Epistatic dominance parameter for alpha_1 in the text
alpha = 0.25
#Strength of dominance.
h =0.5


#Function to calculate probability of fixation

#Helper functions for Kimura's prob of fixation
nu = function(i,s,A){
  if(length(A)>0 ) x = prod(sapply(A, function(x) 1+epsilon(i,x)))
  else x = 1
  return((1+s)*x-1)
}

mu = function(i,s,h,alpha,A){
  x = nu(i,s,A)
  if(length(A)>0) {y = prod(sapply(A, function(x) 1+alpha*epsilon(i,x)))
  } else {y = 1}
  if(x==0){
    return(1/2)
  } else{
    return(((1+h*s)*y-1)/x) 
  }
}
#Calculating the actual probability of fixation
probFix = function(i,A){
  s = Mutations$s[i]
  c = Ne*nu(i,s,A)
  D = 2*mu(i,s,h,alpha,A)-1
  integ = function(x) (exp(-2*c*D*x*(1-x)-2*c*x))
  if(c > -1000){
    x = integrateR(integ ,0,1/(2*Ne),ord=14,verbose=FALSE)$value
    y = integrateR(integ ,0,1,ord=14,verbose=FALSE)$value
    return(x/y)
  } else {
    return(0)
  }
}



#Precompiled version of the function, useful to save time in simulation

probFixC = cmpfun(probFix)

#Calculate Probability of Fixation in the ancestor for Array, Query, and both in two independent populations

probFixA = sapply(data$Array.SMF,function(x) probFix2C(x,1+0.5*(x-1)))
probFixQ = sapply(data$Query.single.mutant.fitness..SMF.,function(x) probFix2C(x,1+0.5*(x-1)))
probFixD = 2*probFixA*probFixQ
probFixD = sapply(data$Array.SMF,function(x) max(x,0))
probFixD2 = sapply(data$Query.single.mutant.fitness..SMF.,function(x) max(x,0))

#Calculate Gamma, which is the probability of observing an interaction in a hybrid times its effect on the hybrid
gamma = probFixD*DmiMult
gamma_Alt = probFixD*data.tot$Genetic.interaction.score..ε.
gamma_Add = probFixD*DmiAdd

#Assign data back

data$ExpFitMult = ExpFitMult
data$ExpFitAdd = ExpFitAdd
data$gamma_Alt = gamma_Alt
data$DMIMult = DmiMult
data$DMIAdd = DmiAdd
data$probFixA = probFixA
data$probFixQ = probFixQ
data$probFix = probFixD
data$gamma = gamma
data$gamma_Add = gamma_Add

#Make sure all is well, overwriting old data now
data.tot = data

#Cleanup of unnecessary vars
rm(data,data.ExE,data.DAmP,data.NxN,data.ExN_NxE)
rm(ExpFitAdd,ExpFitMult,DmiMult,DmiAdd,probFixA,probFixQ,probFixD,gamma,gamma_Add)


#Calculate a summary for a bunch of genes, useful to save time in later analyses, not optimized itself, so takes a while.

StrainData = data.frame(nrow = uniStrains,ncol=7)
for(i in 1:nStrains){
  StrainData[i,1] = uniStrains[i]
  fitness = strainFits[i]
  pF = max(c(probFixC(fitness),0),na.rm=TRUE)
  id1 = which(SFit$Strain.ID == uniStrains[i])
  StrainData[i,2] = SFit$Systematic.gene.name[id1]
  idx1 = which(data.tot$Query.Strain.ID == uniStrains[i])
  idx2 = which(data.tot$Array.Strain.ID == uniStrains[i])
  GI = mean(c(data.tot$DMIMult[idx1],data.tot$DMIMult[idx2]),na.rm=TRUE)
  GI_Add = mean(c(data.tot$DMIAdd[idx1],data.tot$DMIAdd[idx2]),na.rm=TRUE)
  GI_C = mean(c(data.tot$Genetic.interaction.score..ε.[idx1],data.tot$Genetic.interaction.score..ε.[idx2]),na.rm=TRUE)
  StrainData[i,3] = pF
  StrainData[i,4] = GI
  StrainData[i,5] = mean(c(data.tot$gamma[idx1],data.tot$gamma[idx2]),na.rm=TRUE)
  StrainData[i,6] = fitness
  StrainData[i,7] = GI_Add
  StrainData[i,8] = GI_C
}
names(StrainData) = c("Strain","Gene","ProbFix","Epsilon","Gamma","Fitness","Epsilon_Add","Epsilon_C")

rm(fitness,pF,id1,idx1,idx2,GI,GI_Add)

StrainData$Epsilon_hi = StrainData$s

i1 = which(data.tot$Query.Strain.ID %in% StrainData$id[highest_obs])
i2 = which(data.tot$Array.Strain.ID %in% StrainData$id[highest_obs])
highest_int = data.tot[unique(c(i1,i2)),]

for(i in 1:nStrains){
  idx1 = which(data.tot$Query.Strain.ID == StrainData$id[i])
  idx2 = which(data.tot$Array.Strain.ID == StrainData$id[i])
  StrainData$Epsilon[i] = mean(c(data.tot$DMIMult[idx1],data.tot$DMIMult[idx2]),na.rm=TRUE)
  idx1 = which(highest_int$Query.Strain.ID == StrainData$id[i])
  idx2 = which(highest_int$Array.Strain.ID == StrainData$id[i])
  StrainData$Epsilon_hi[i] = mean(c(highest_int$DMIMult[idx1],highest_int$DMIMult[idx2]),na.rm=TRUE)
}
  
### Simulation ###


#Next start on Simulation. We first make everything faster by precomputing a lookup table for epistatic effect strength, called EpiTable
#Here's the table for simulations using the epistasis metric used in the Costanzo paper. Change out the epsilon parameters being added for other epistasis.
EpsiTable = matrix(data = NA,nrow = nStrains,ncol = nStrains)

for(i in 1:nStrains){
  idxA = which(data.tot$Query.Strain.ID == uniStrains[i])
  idxB = which(data.tot$Array.Strain.ID == uniStrains[i])
  for(j in i:nStrains){
    idxA2 = which(data.tot$Array.Strain.ID[idxA] == uniStrains[j])
    idxB2 = which(data.tot$Query.Strain.ID[idxB] == uniStrains[j])
    epi1 = data.tot$Genetic.interaction.score..ε.[idxA][idxA2]
    epi2 = data.tot$Genetic.interaction.score..ε.[idxB][idxB2]
    l = c(0,abs(epi1),abs(epi2))
    idx = which.max(l)
    EpsiTable[i,j] = c(0,epi1,epi2)[idx]
  }
}
save("epiTableCostanzo.RData",EpsiTable)
rm(idxA2,idxB2,epi1,epi2,l,idx)



load("epiTableCostanzo.RData")

#Define a function which will randomly sample which allele is fixed, given a list of probabilities:


probFix2 = function(nu,mu){
  if(nu<0.25 & mu < 0.25){
    return(0)
  }
  c = Ne*(nu-1)
  if(c!= 0) {
    D = 2*(mu-1)/(nu-1)-1
  } else {
    D = 0
  }
  integ = function(x) (exp(-2*c*D*x*(1-x)-2*c*x))
  if(c > -100 & c*D > -10000){
    x = integrateR(integ ,0,1/(2*Ne),ord=4,verbose=FALSE)$value
    y = integrateR(integ ,1/(2*Ne),1,ord=4,verbose=FALSE)$value
    return(max(x/(x+y),0))
  } else {
    return(0)
  }
}

probFix2C = compile(probFix2)

nextFix2= function(probs){
  norm = probs/sum(probs,na.rm=TRUE)
  norm[is.na(norm)] = 0
  l = length(probs)
  i = sample(1:l,1,prob=norm)
  return(i)
}

nextFix2C = compile(nextFix2)

#The simulation: runs for K substitutions, given alpha1, alpha2 and h
fixSimu2 = function(K,alpha1,alpha2,h){
  alpha = alpha2 #This is just to simplify the dependence between parameters in debugging
  nuA = Mutations$s+1
  muA = rep(h,nStrains)*Mutations$s+1
  canFixA = 1:nStrains
  canFixB = 1:nStrains
  nuB = Mutations$s+1
  muB = rep(h,nStrains)*Mutations$s+1
  fA = c()
  fB = c()
  for(i in 1:K){
    #Alternately fix a substitution in A, then in B
    #PopA
    probsA = sapply(canFixA,function(x) probFix2C(nuA[x],muA[x]))
    if(sum(probsA,na.rm=TRUE)>0){
      nA = canFixA[nextFix2(probsA)]
      fA = c(fA,nA)
      canFixA = setdiff(canFixA,nA)
      canFixB = setdiff(canFixB,nA)
      nuA = nuA*(1+c(EpsiTable[(1:nA),nA],EpsiTable[nA,(nA+1):nStrains]))
      muA = muA*(1+alpha*c(EpsiTable[(1:nA),nA],EpsiTable[nA,(nA+1):nStrains]))
      muA[which(muA<0)] = 0
      nuA[which(nuA<0)] = 0
    }  else {
      fA = c(fA,rep(NA,(K-i)))
      fB = c(fB,rep(NA,(K-i)))
      return(data.frame(fA=fA,fB=fB,check.rows=FALSE)) 
    }
    #PopB
    probsB = sapply(canFixB,function(x) probFix2C(nuB[x],muB[x]))
    if(sum(probsB,na.rm=TRUE)>0){
      nB = canFixB[nextFix2(probsB)]
      fB = c(fB,nB)
      canFixA = setdiff(canFixA,nB)
      canFixB = setdiff(canFixB,nB)
      nuB = nuB*(1+c(EpsiTable[(1:nB),nB],EpsiTable[nB,(nB+1):nStrains]))
      muB = muB*(1+alpha*c(EpsiTable[(1:nB),nB],EpsiTable[nB,(nB+1):nStrains]))
    } else {
      fA = c(fA,rep(NA,(K-i)))
      fB = c(fB,rep(NA,(K-i)))
      return(data.frame(fA=fA,fB=fB,check.rows=FALSE)) 
    }
  }
  return(data.frame(fA=fA,fB=fB,check.rows = FALSE))
}

fixSimu2C = compile(fixSimu2)


#The Actual Simulation, takes the total number of mutations to fix (K) as a parameter, and fixes mutations in alternating populations until K are fixed.

#Running the Simulations
#I have made some minimal parallelization effort to these. It's not necessary - to avoid just use a %do% foreach loop instead of %dopar%


#Set pop sizes, effective and actual
Ne = 10e06
N = 10e06
#Number of substitutions per population
K = 50


cl = makeCluster(6) #Number of cores. I find that with 16 cores I can run around 8 effectively, since each cluster will multithread further. 6 for when I want to be able to use laptop.
registerDoSNOW(cl) #Initialize cluster
start = 1 #Used to continue runs I cancel for whatever reason.
iterations = 10000 #Number of runs desired. End goal is 100,00, else wouldn't bother with multithreading.
done = sapply(list.files(path="./SimuRes2018/"),function(x) as.numeric(sub(".csv","",sub("ret_","",x))))
RES = foreach(i=setdiff(start:iterations,done),.combine=cbind,.verbose=T) %dopar% { #Will return a matrix of the results combined with cbind. Verbose helps keep track/debug, not necessary.
  require(Rmpfr)
  #Change the parameters here for alternate dominance scenarios, but remember to use the same parameters in the analyses below, they are not saved with the simulation result
  ret = fixSimu2C(K,0.25,0.5,0.5) #Actually run the simulation
  write.csv(ret,file=paste("./SimuRes2018/ret_",i,".csv",sep="")) #Save the results in a file. Causes some slowdown, but simulations take long enough you should be saving these anyways.
  return(ret)
}

stopCluster(cl) #Stop cluster

##Reading simulation results back in.

FA_2 = matrix(data=NA, nrow=1000,ncol=50)
FB_2 = matrix(data=NA,nrow=1000,ncol=50)
for(i in 1:1000){
  ret = read.csv(paste("./SimuRes2018/ret_",i,".csv",sep=""),header=TRUE,stringsAsFactors = FALSE)
  FA_2[i,] = ret$fA
  FB_2[i,] = ret$fB
}

# Reading data back in.

FA = matrix(data=NA,nrow=iterations,ncol=ceiling(K))
FB = matrix(data=NA,nrow=iterations,ncol=ceiling(K))
for(i in 1:iterations){
  ret = read.csv(paste("./SimuRes2018/ret_",i,".csv",sep=""),header=TRUE,stringsAsFactors = FALSE)
  FA[i,] = ret$fA
  FB[i,] = ret$fB
}

#Calculating frequency of strains:
test = table(c(c(FA),c(FB)))/iterations
idx  = names(test)
test2 = c()
for(i in 1:nStrains){
  n = test[which(names(test)==i)]
  if(length(n)==0){test2 = c(test2,0)} else {test2 = c(test2,as.numeric(n))}
}
StrainData$ObsFix = test2


rm(test,test2)


###Analysis
#Define several functions that calculate parental and hybrid fitness, etc.
#First define the dominance of epistasis in a heterozygote for both loci
alpha2 = 0.5
#calculates fitness of derived population with set fA of fixed alleles.
calcFit = function(fA,alpha1=0.25){
  dir = prod(1+StrainData$s[fA])
  #Epistatic interaction effects
  epi = prod(1+EpsiTable[fA,fA],na.rm=TRUE)
  return(max(dir*epi,0))
}

#Calculate Relative Hybrid Fitness, given fixed mutations in pop A (fA), and B (fB)

calcHybrFit = function(fA,fB,alpha=0.5,h=0.5){
  alpha1 = alpha2*alpha
  #All interactions (ancestral and new) are at 1/4 strength, calculate those
  epi = prod(1+alpha1*EpsiTable[c(fA,fB),c(fA,fB)],na.rm=TRUE)
  #All direct fitness effects are at heterozygous state
  dir = prod(1+h*(StrainData$s[c(fA,fB)]),na.rm=TRUE)
  #Return as a function of the fitness of parental populations
  return(epi*dir/(0.5*calcFit(fA)+0.5*calcFit(fB)))
}
calcHybrFit_rec = function(fA,fB,alpha=0,h=0.5){
  alpha1 = alpha2*alpha
  #All interactions (ancestral and new) are at 1/4 strength, calculate those
  epi = prod(1+alpha1*EpsiTable[c(fA,fB),c(fA,fB)],na.rm=TRUE)
  #All direct fitness effects are at heterozygous state
  dir = prod(1+h*(StrainData$s[c(fA,fB)]),na.rm=TRUE)
  #Return as a function of the fitness of parental populations
  return(epi*dir/(0.5*calcFit(fA)+0.5*calcFit(fB)))
}
calcHybrFit_dom = function(fA,fB,alpha=1,h=0.5){
  alpha1 = alpha2*alpha
  #All interactions (ancestral and new) are at 1/4 strength, calculate those
  epi = prod(1+alpha1*EpsiTable[c(fA,fB),c(fA,fB)],na.rm=TRUE)
  #All direct fitness effects are at heterozygous state
  dir = prod(1+h*(StrainData$s[c(fA,fB)]),na.rm=TRUE)
  #Return as a function of the fitness of parental populations
  return(epi*dir/(0.5*calcFit(fA)+0.5*calcFit(fB)))
}
calcHybrFitMax = function(fA,fB,alpha1=0.25,h=0.5){
  #All interactions (ancestral and new) are at 1/4 strength, calculate those
  epi = prod(1+alpha1*EpsiTable[c(fA,fB),c(fA,fB)],na.rm=TRUE)
  #All direct fitness effects are at heterozygous state
  dir = prod(1+h*(StrainData$s[c(fA,fB)]),na.rm=TRUE)
  #Return as a function of the fitness of parental populations
  return(epi*dir/max(c(calcFit(fA),calcFit(fB))))
}
calcHybrFitMin = function(fA,fB,alpha1=0.25,h=0.5){
  #All interactions (ancestral and new) are at 1/4 strength, calculate those
  epi = prod(1+alpha1*EpsiTable[c(fA,fB),c(fA,fB)],na.rm=TRUE)
  #All direct fitness effects are at heterozygous state
  dir = prod(1+h*(StrainData$s[c(fA,fB)]),na.rm=TRUE)
  #Return as a function of the fitness of parental populations
  return(epi*dir/min(c(calcFit(fA),calcFit(fB))))
}
#Calculating the average epistatic effect within a population
calcEps = function(fA){
  return(mean(EpsiTable[fA,fA],na.rm=TRUE))
}
#Calculating the average epistatic effect between populations
calcEpsHybr = function(fA,fB){
  c1 = c(EpsiTable[fA,fB],EpsiTable[fB,fA])
  return(mean(c1,na.rm=TRUE))
}
#Calculate the NUMBER of interactions in a hybrid
calcNumEps = function(fA,fB){
  with = length(intersect(which(EpsiTable[fA,fA]!=0),which(!is.na(EpsiTable[fA,fA]))))+
    length(intersect(which(EpsiTable[fB,fB]!=0),which(!is.na(EpsiTable[fB,fB]))))
  betw = length(intersect(which(EpsiTable[fA,fB]!=0),which(!is.na(EpsiTable[fA,fB]))))+
    length(intersect(which(EpsiTable[fB,fA]!=0),which(!is.na(EpsiTable[fB,fA]))))
  return(c(with,betw))
}
calcNumEpsBetw = function(fA,fB){
  betw = length(intersect(which(EpsiTable[fA,fB]!=0),which(!is.na(EpsiTable[fA,fB]))))+
    length(intersect(which(EpsiTable[fB,fA]!=0),which(!is.na(EpsiTable[fB,fA]))))
  return(betw)
}
calcNumEpsWith = function(fA,fB){
  with = length(intersect(which(EpsiTable[fA,fA]!=0),which(!is.na(EpsiTable[fA,fA]))))+
    length(intersect(which(EpsiTable[fB,fB]!=0),which(!is.na(EpsiTable[fB,fB]))))
  return(with)
}

#Calculate the effective selection coefficient for mutation i in the background of fB fixed mutations
calcSe = function(i,fB){
  epi = prod(1+EpsiTable[i,fB],na.rm=TRUE)*prod(1+EpsiTable[fB,i])
  return((1+StrainData$s[i])*epi)
}
#Calculate the fraction of substitutions in fA that are positive in the background of fB
calcFracPos = function(fA,fB){
  retA = c()
  retB = c()
  for(i in 1:length(fA)){
    retA = c(retA,length(which(sapply(1:i,function(x) calcSe(fA[x],fB[1:i]))>1))/i)
    retB = c(retB,length(which(sapply(1:i,function(x) calcSe(fB[x],fA[1:i]))>1))/i)
  }
  return(sapply(1:length(fA),function(x) mean(c(retA[x],retB[x]))))
}
#Calculate the ratio of substitutions from A that can introgress to B, vs substitutions from B that can introgress to A
calcFracPosDir =  function(fA,fB){
  retA = c()
  retB = c()
  for(i in 1:length(fA)){
    retA = c(retA,length(which(sapply(1:i,function(x) calcSe(fA[x],fB[1:i]))>1))/i)
    retB = c(retB,length(which(sapply(1:i,function(x) calcSe(fB[x],fA[1:i]))>1))/i)
  }
  return(sapply(1:length(fA),function(x) mean(retA[x])/mean(retB[x])))
}
#Helper functions to loop over the sets and calculate fitness/epistasis over varying timespans
doForSeq = function(v,f)  {
  #Takes vector v and applies function f to sequentially larger subsets of v.
  m = length(v)
  ret = c()
  for(i in 1:m){
    ret = c(ret,f(v[1:i]))
  }
  return(ret)
}
doFor2Seq = function(v1,v2,f)  {
  #same as do for seq, but for functions with 2 parameters, and takes two sets of equal size
  m = length(v1)
  ret = c()
  for(i in 1:m){
    ret = c(ret,f(v1[1:i],v2[1:i]))
  }
  return(ret)
}

#Performing the analyses. This takes a while
#First define your dominance alpha_2
alpha2 = 0.95
fit_resA = apply(FA,1,function(x) doForSeq(x,calcFit))
fit_resB = apply(FB,1,function(x) doForSeq(x,calcFit))
epi_resA = apply(FA,1,function(x) doForSeq(x,calcEps))
epi_resB = apply(FB,1,function(x) doForSeq(x,calcEps))

frac_pos = sapply(1:iterations,function(x) calcFracPos(FA[x,],FB[x,]))
frac_pos_dir = sapply(1:iterations,function(x) calcFracPosDir(FA[x,],FB[x,]))
epi_resHyb = sapply(1:iterations,function(x) doFor2Seq(FA[x,],FB[x,],calcEpsHybr))
fit_resHyb = sapply(1:iterations,function(x) doFor2Seq(FA[x,],FB[x,],calcHybrFit))
fit_resHyb_dom = sapply(1:iterations,function(x) doFor2Seq(FA[x,],FB[x,],calcHybrFit_dom))
fit_resHyb_rec = sapply(1:iterations,function(x) doFor2Seq(FA[x,],FB[x,],calcHybrFit_rec))
fit_resHybMax = sapply(1:iterations,function(x) doFor2Seq(FA[x,],FB[x,],calcHybrFitMax))
fit_resHybMin = sapply(1:iterations,function(x) doFor2Seq(FA[x,],FB[x,],calcHybrFitMin))

num_epiHybWith = sapply(1:iterations,function(x) doFor2Seq(FA[x,],FB[x,],calcNumEpsWith))
num_epiHybBetw = sapply(1:iterations,function(x) doFor2Seq(FA[x,],FB[x,],calcNumEpsBetw))

fit_within = cbind(fit_resA,fit_resB)
epi_within = cbind(epi_resA,epi_resB)

#Calculating frequency of interactions being observed within or between populations (cooc vs sep). Must be a more efficient way.
#

coocTables = matrix(nrow =nStrains,ncol=nStrains)
sepTables = matrix(nrow=nStrains,ncol=nStrains)
for(i in 1:nStrains){
  idxA = apply(FA,1,function(x) which(x == i))
  idxA = as.numeric(idxA[sapply(idxA,length)>0])
  idxB = apply(FB,1,function(x) which(x == i))
  idxB = as.numeric(idxB[sapply(idxB,length)>0])
  for(j in i:nStrains){
    idxA2 = length(which(FA[idxA,] == j))
    idxB2 = length(which(FB[idxB,] == j))
    coocTables[i,j] = idxA2+idxB2
    idxA2 = length(which(FB[idxA,] == j))
    idxB2 = length(which(FA[idxB,] == j))
    sepTables[i,j] = idxA2+idxB2
  }
}
sepTables = sepTables/iterations
coocTables = coocTables/iterations

save(sepTables,file="SepTables.Rdata")
save(coocTables,file="CoocTables.Rdata")

#annotating data.tot with cooc and sep values
for(i in 1:nStrains){
  idx1 = which(data.tot$Query.Strain.ID == Mutations$id[i])
  idx2 = match(data.tot$Array.Strain.ID[idx1],Mutations$id)
  idx1.2 = idx1[which(idx2>i)]
  idx2.2 = idx2[which(idx2>i)]
  data.tot$Cooc.Obs[idx1.2] = coocTables[i,idx2.2]
  data.tot$Sep.Obs[idx1.2] = sepTables[i,idx2.2]
  idx1 = which(data.tot$Array.Strain.ID == Mutations$id[i])
  idx2 = match(data.tot$Query.Strain.ID[idx1],Mutations$id)
  idx1.2 = idx1[which(idx2>i)]
  idx2.2 = idx2[which(idx2>i)]
  data.tot$Cooc.Obs[idx1.2] = coocTables[i,idx2.2]
  data.tot$Sep.Obs[idx1.2] = sepTables[i,idx2.2]
}

#Calculating the gamma_obs value (gamma as a function of observed frequency of interactions instead of expected.)
# formula is equal to: gamma_n(i,j)  = epsilon_ij*{(Frequncy Observed between/n^2)-(Frequency observed within/(n(n-1)))}
n = ceiling(K/2)
gammaObs = data.tot$DMIMult*(data.tot$Sep.Obs/(n^2)-data.tot$Cooc.Obs/(n*(-1)))

data.tot$gamma_50 = gammaObs                     
                             
#PLOTS
library(ggplot2)
library(hexbin)
library(cowplot)



#Figure: a)x-> Fitness Strain 1, y-> fitness strain2, color-> average epistatic effect 
#b) Distribution of direct fitness effects. 
#c) Distribution of epistatic effects. 
#d) Zoom in of region in a with likely to evolve genes (fitness > 1) 
#e) Relationship between joint likelihood of fixation (x) and the average strength of DMI (y)

par(cex = 1.5)
layout(matrix(c(1,2,3,4,4,4),2,3,byrow=TRUE))
#1a)
dens = density(StrainData$Fitness,na.rm=TRUE)
plot(dens,main=NA,xlab=NA,ylab=NA,lwd=2)
x1 = min(which(dens$x >= 0))
x2 = max(which(dens$x <= 1.2))
with(dens, polygon(x=c(x[c(x1:x2)]),y=c(y[x1:x2]),col=rgb(0,0,0,0.2)))

plot1a = ggplot(StrainData,aes(x=s)) + 
  geom_density(fill="grey") + 
  labs(x="Selection Coefficient",y="Density") + 
  geom_vline(xintercept=median(StrainData$s,na.rm=TRUE),color="red")+
  scale_x_continuous(breaks=seq(-1,0,0.25))


#1b)
dens = density(data.tot$DMIMult/4,na.rm=TRUE)
x1 = min(which(dens$x >= -0.3))
x2 = max(which(dens$x <= 0.2))
plot(dens,main=NA,xlab=NA,ylab=NA,xlim=c(-0.3,0.3))
with(dens, polygon(x=c(x[c(x1:x2)]),y=c(y[x1:x2]),col=rgb(0,0,0,0.2)))
abline(v=0,col="red",lwd = 3)

plot1b = ggplot(data.tot,aes(x=DMIMult)) + 
  geom_density(fill="grey") + 
  geom_vline(xintercept = median(data.tot$DMIMult,na.rm=TRUE),color="red") + 
  coord_cartesian(xlim=c(-0.3,0.3)) + 
  labs(x="Epistatic Effect",y="Density")


#c) Not the best way to do this plot, but it does show one of the main points - most of the space is close to 0 (i.e. white)
#Average Epistatic effect vs. Fitness
plot(StrainData$Fitness,StrainData$Epsilon/4,pch=20,col=rgb(0,0,0,0.1),xlab=NA,ylab=NA)
abline(h=0,col="red",lwd=2)

plot1c = ggplot(StrainData,aes(x=s,y=Epsilon)) + 
  geom_point(alpha = 0.1) + 
  geom_hline(yintercept = 0,color="black",linetype=2)  + 
  geom_smooth(method="lm",aes(color="red"),show.legend=FALSE) + 
  coord_cartesian(ylim=(c(0.05,-0.05))) + 
  labs(x="Selection Coefficient",y="Mean Epistasis")+
  scale_x_continuous(breaks=seq(-1,0,0.25))


asinh_trans <- function(){
  trans_new(name = 'asinh', transform = function(x) asinh(x), 
            inverse = function(x) sinh(x))
}

plot1d = ggplot(data.tot,aes(x=Query.single.mutant.fitness..SMF.-1,y= Array.SMF-1,z=DMIMult)) + 
  stat_summary_2d(bins=30,fun ="mean",show.legend = TRUE) + 
  scale_fill_gradient2(name="Mean",high="red",low="blue",mid="gray96",midpoint = 0) + 
  labs(x="Selection Coefficient 1",y="Selection Coefficient 2") +
  geom_vline(xintercept=0,color="red")+geom_hline(yintercept=0,color="red")


abline(h=1,v=1,lty=2)

#d) Same problem as a, replace stat_summary_hex with geom_point for better visual
idx1 = which(data.tot$Query.single.mutant.fitness..SMF.>=1) 
idx2 = which(data.tot$Array.SMF >=1)
sigData = data.tot[intersect(idx1,idx2),]
plot1e = ggplot(sigData,aes(x=Query.single.mutant.fitness..SMF.-1 ,y= Array.SMF-1,z=DMIMult)) +
  stat_summary_2d(bins=20,fun = "median",show.legend = TRUE) + 
  scale_fill_gradient2(high="red",low="blue",mid="gray96",midpoint = 0) + 
  labs(x="Direct Fitness Effect 1",y="Direct Fitness Effect 2")


percent_Pos = function(x){
  idx = which(data.tot$P.value<x)
  idx2 = which(data.tot[idx,]$DMIMult<0)
  idx3 = which(data.tot[idx,]$DMIMult>0)
  return(length(idx3)/(length(idx2)+length(idx3)))
}
#with epsilon cutoff
percent_Pos = function(p,epsilon){
  idx = intersect(which(data.tot$P.value<p),which(abs(data.tot$DMIMult)>epsilon))
  idx2 = which(data.tot[idx,]$DMIMult<0)
  idx3 = which(data.tot[idx,]$DMIMult>0)
  return(length(idx3)/(length(idx2)+length(idx3)))
}

#Figure to show proportion of positive interactions
posIntProp = sapply(seq(0.001,0.2,0.005),function(x) percent_Pos(0.05,x))

plot1b_new = ggplot()+
  geom_line(aes(x=seq(0.001,0.2,0.005),y=posIntProp),col="red",size=2)+
  coord_cartesian(ylim=c(0,0.5))+
  labs(x="Interaction Strength Cutoff",y="Proportion Positive Interactions")

#e) We remove the gamma ~= 0 points to reduce number of drawn points drastically
idx1 = intersect(which(data.tot$P.value<0.5),which(data.tot$ProbFix > 0))
plot1f = ggplot(data.tot[idx1,],aes(x=ProbFix,y=DMIMult,z=DMIMult))+
  #geom_rug(show.legend=FALSE,aes(color="black"))+
  stat_summary_hex(fun="length",show.legend = TRUE)+
  scale_fill_gradientn(name="Log(Occurrence)",trans="log",colours=rainbow(10),breaks=c(1,10,100,1000,10000,100000))+
  labs(x="Probability of Interaction",y="Epistatic effect")



#Densities if needed
plot1fx = ggplot(data.tot[idx1,],aes(x=probFix))+geom_density(fill="grey")
plot1fy = ggplot(data.tot[idx1,],aes(x=DMIMult))+geom_density(fill="grey")


#fig1 = plot_grid(plot1a,plot1b,plot1c,plot1d,plot1e,NULL,plot1f,NULL,NULL,nrow=3,align="h",labels = c("A","B","C","D","E","F"))
fig1 = ggdraw()+draw_plot(plot1b,0, 3/4, 1/2, 1/4) +
 draw_plot(plot1b_new,1/2, 3/4, 1/2, 1/4) +
 draw_plot(plot1c,0, 1/2, 1/2, 1/4) +
 draw_plot(plot1d,1/2, 1/2, 1/2, 1/4) +
 draw_plot(plot1f,1/5, 0, 3/4, 1/2) +
 draw_plot_label(c("A","B","C","D","E"),c(0,1/2,0,1/2,1/5),c(1,1,3/4,3/4,1/2),size=14)


fig2 = ggdraw()+draw_plot(plot1b_new,0,3/4,1/2,1/4)+
  draw_plot(plot1c,1/2,3/4,1/2,1/4)+
  draw_plot(plot1f,0,0,1,3/4)

                                                                                                                                                                                                
#Figure 2: a)Expected vs. Observed Fixation Rates b) Hybrid Fitness over time c) Between vs within pop epsilon

#2a

par(cex = 1.5)
plot(StrainData$Fitness-1,StrainData$ObsFix50,xlab="Fitness Coefficient of Mutation",ylab="Observed Fixation Rate")


plot2a = ggplot(StrainData,aes(x=s,y=ObsFix,color=Epsilon)) + 
  geom_point() + scale_colour_gradient2(high="red",low="green",mid="yellow",midpoint=0)

plot2a = ggplot(StrainData) + geom_point(aes(x=s,y=ObsFix,color=Epsilon))+
  scale_color_gradient2(high="red",mid="grey90",low="blue",midpoint=0)+
  labs(x="Selection Coefficient",y="Observed Fixation Rate") 

plot2a =  ggplot() + 
  geom_point(data = subset(StrainData,Epsilon < 0.0),aes(x=s,y=ObsFix),color=wes_palette("Darjeeling1")[1],alpha=0.25) +
  geom_point(data = subset(StrainData,Epsilon > 0.0),aes(x=s,y=ObsFix),color=wes_palette("Darjeeling1")[2],alpha=0.25) + 
  labs(x="Selection Coefficient",y="Observed Fixation Rate") 


#2b - non ggplot way
par(cex=1.5)
matplot(fit_resHyb_rec[,100:200],type="l",lty=1,col=rgb(0,0,0,0.1),ylim=c(0,1.25),xlab = "Fixed Substitutions",ylab="Relative Hybrid Fitness")
lines(apply(fit_resHyb_rec,1,mean,na.rm=TRUE),lwd=3,lty=2,col="red",type="l",ylim=c(0,1.25),xlab="Fixed Substitutions",ylab="Relative Hybrid Fitness")
coords.x = c(1:K,rev(1:K))
upper = apply(fit_resHyb_rec,1,quantile,na.rm=TRUE,probs=0.975)
lower = apply(fit_resHyb_rec,1,quantile,na.rm=TRUE,probs=0.025)
coords.y = c(upper,rev(lower))
polygon(coords.x,coords.y,col=rgb(0.25,0.25,0.25,0.25),border=NA)
abline(h=1,col="black",lwd=3)

fitDF = melt(fit_resHyb_rec[,500:600])
fitDF_dom = melt(fit_resHyb_dom[,500:600])
fitDF_rec = melt(fit_resHyb_rec[,500:600])
l1 = exp(1/4*(en+eo)*(1:K)-choose(1:K,2)*eo*1/2)
l2 = exp(et*(choose(1:K,2)))

fitDFAdd = fitDF
fitDFAdd$values = as.vector(fit_resHybAdd[,1:100])
fitDFC = fitDF
fitDFC$values = as.vector(fit_resHybC[,1:100])

#ggplot way
plot2b = ggplot() + geom_ribbon(aes(x=1:50,ymin=lower,ymax=upper),alpha=0.5,fill="darkgoldenrod3",col="darkgoldenrod4") + 
  labs(x="Fixed Substitutions (per population)",y="Relative Hybrid Fitness") + 
  geom_line(aes(x=1:50,y=apply(fit_resHyb_rec,1,mean,na.rm=TRUE)),size=1.5)+
  geom_hline(yintercept = 1,color="grey50",linetype="longdash",size=0.5) + 
  #geom_line(aes(x=1:50,y=l1),size=2,color="blue") + 
  geom_line(data = fitDF,aes(x=Var1,y=value,group=Var2),alpha=0.25) + 
  #geom_line(aes(x=1:50,y=l2),size=2,color="purple") +
  coord_cartesian(xlim=c(1,50),ylim=c(0,1.1),expand=FALSE)


plot2b


plots2b = ggplot() +
  labs(x="Fixed Substitutions (per population)",y="Relative Hybrid Fitness") + 
  geom_hline(yintercept = 1,color="grey50",linetype="longdash",size=0.5) + 
  #geom_line(aes(x=1:50,y=l1),size=2,color="blue") + 
  #geom_line(aes(x=1:50,y=l2),size=2,color="purple") +
  coord_cartesian(xlim=c(1,50),ylim=c(0,2.0),expand=FALSE)+
  geom_ribbon(aes(x=1:50,
                  ymin=apply(fit_resHyb,1,quantile,na.rm=TRUE,probs=0.025),
                  ymax=apply(fit_resHyb,1,quantile,na.rm=TRUE,probs=0.975)),
              alpha=0.25,fill="darkgoldenrod3",col="darkgoldenrod4") +
  
  geom_ribbon(aes(x=1:50,
                  ymin=apply(fit_resHyb_rec,1,quantile,na.rm=TRUE,probs=0.025),
                  ymax=apply(fit_resHyb_rec,1,quantile,na.rm=TRUE,probs=0.975)),
              alpha=0.25,fill="slateblue1",col="slateblue2") +
  
  geom_ribbon(aes(x=1:50,
                  ymin=apply(fit_resHyb_dom,1,quantile,na.rm=TRUE,probs=0.025),
                  ymax=apply(fit_resHyb_dom,1,quantile,na.rm=TRUE,probs=0.975)),
              alpha=0.25,fill="tomato2",col="tomato3") +
  geom_line(aes(x=1:50,y=apply(fit_resHyb,1,mean,na.rm=TRUE)),size=1.5,col="darkgoldenrod3")+
  geom_line(aes(x=1:50,y=apply(fit_resHyb_rec,1,mean,na.rm=TRUE)),size=1.5,col="slateblue1")+
  geom_line(aes(x=1:50,y=apply(fit_resHyb_dom,1,mean,na.rm=TRUE)),size=1.5,col="tomato2")

#2c) non-ggplot, then ggplot
plot(1:K,rep(0,K),lty = 3,col=rgb(0,0,0,0.75),lwd = 5,xlim=c(4,50),ylim=c(-0.002,0.002),type="l",xlab="Fixed Substitutions",ylab="Average Epistatic Effect")

plot2c = ggplot()+
  geom_ribbon(aes(x=2:K,ymin=lower1[2:K],ymax=upper1[2:K]),fill="darkslateblue",alpha=0.5)+
  geom_ribbon(aes(x=2:K,ymin=lower2[2:K],ymax=upper2[2:K]),fill="darkgoldenrod",alpha=0.5) + 
  coord_cartesian(xlim=c(2,50),expand=FALSE) + 
  geom_hline(yintercept = 0,linetype = "longdash",size=0.5,col="grey50") + 
  labs(x="Fixed Substitutions",y="Average Epistatic Effect") + 
  geom_line(aes(x=2:K,y=apply(epi_resHyb[2:50,],1,mean,na.rm=TRUE)),size=1.5,color="darkgoldenrod") + 
  geom_line(aes(x=2:K,y=apply(epi_within[2:50,],1,mean,na.rm=TRUE)),size=1.5,color="darkslateblue")

#Want to plot +-stdev, first for within, then across pop epsilon, find it easier to precalculate these values
coords.x = c(1:K,rev(1:K))

upper1 = apply(epi_within,1,quantile,na.rm=TRUE,probs=0.975)
lower1 = apply(epi_within,1,quantile,na.rm=TRUE,probs=0.025)
coords.y = c(upper1,rev(lower1))
polygon(coords.x,coords.y,col=rgb(0,0,1,0.25),border=NA)


upper2 = apply(epi_resHyb,1,quantile,na.rm=TRUE,probs=0.975)
lower2 = apply(epi_resHyb,1,quantile,na.rm=TRUE,probs=0.025)
coords.y = c(upper2,rev(lower2))
polygon(coords.x,coords.y,col=rgb(1,0,0,0.25),border=NA)

lines(apply(epi_resHyb,1,mean,na.rm=TRUE),col=rgb(1,0,0,1),lwd=2,lty=1)
lines(apply(epi_within,1,mean,na.rm=TRUE),col=rgb(0,0,1,1),lwd=2,lty=1)

fig2 = ggdraw()+draw_plot(plot2a,0, 2/3, 1/2, 1/3) +
  draw_plot(plot2b,0, 0, 1, 2/3) +
  draw_plot(plot2c,1/2, 2/3, 1/2, 1/3) +
  draw_plot_label(c("A","C","B"),c(0,0,1/2),c(1,2/3,1),size=15)

plot(epi_within[,1],ylim = c(-0.01,0.02))
for(i in 2:1000){
  lines(epi_within[,i],lwd=0.25,col="grey")
}


#Figure 3: Change in the importance of interactions over time.

coords.x = c(1:K,rev(1:K))

upper1 = apply(num_epiHybWith,1,quantile,na.rm=TRUE,probs=0.95)
lower1 = apply(num_epiHybWith,1,quantile,na.rm=TRUE,probs=0.05)
upper2 = apply(num_epiHybBetw,1,quantile,na.rm=TRUE,probs=0.95)
lower2 = apply(num_epiHybBetw,1,quantile,na.rm=TRUE,probs=0.05)
plot3 = ggplot()+
  geom_ribbon(aes(x=2:K,ymin=lower1[2:K],ymax=upper1[2:K]),fill=wes_palette("Darjeeling1")[2],alpha=0.5)+
  geom_ribbon(aes(x=2:K,ymin=lower2[2:K],ymax=upper2[2:K]),fill=wes_palette("Darjeeling1")[1],alpha=0.5) + 
  coord_cartesian(xlim=c(2,50),expand=FALSE) + 
  labs(x="Fixed Substitutions",y="Number of Interactions") + 
  geom_line(aes(x=2:K,y=apply(num_epiHybBetw[2:50,],1,mean,na.rm=TRUE)),size=1.5,color=wes_palette("Darjeeling1")[1]) + 
  geom_line(aes(x=2:K,y=apply(num_epiHybWith[2:50,],1,mean,na.rm=TRUE)),size=1.5,color=wes_palette("Darjeeling1")[2])

#The other way to look at this is as a fraction of possible interactions, this is more interesting
expected = sapply(1:50,function(x) x*(x-1))*2
expected2 = sapply(1:50,function(x) x*(x))

upper1 = apply(num_epiHybWith,1,quantile,na.rm=TRUE,probs=0.95)/expected
lower1 = apply(num_epiHybWith,1,quantile,na.rm=TRUE,probs=0.05)/expected
upper2 = apply(num_epiHybBetw,1,quantile,na.rm=TRUE,probs=0.95)/expected2
lower2 = apply(num_epiHybBetw,1,quantile,na.rm=TRUE,probs=0.05)/expected2
plot3 = ggplot()+
  geom_ribbon(aes(x=2:K,ymin=lower1[2:K],ymax=upper1[2:K]),fill=wes_palette("Darjeeling1")[1],alpha=0.5)+
  geom_ribbon(aes(x=2:K,ymin=lower2[2:K],ymax=upper2[2:K]),fill=wes_palette("Darjeeling1")[2],alpha=0.5) + 
  coord_cartesian(xlim=c(2,50),expand=FALSE) + 
  labs(x="Fixed Substitutions",y="Proportion of possible interactions") + 
  geom_line(aes(x=2:K,y=sapply(2:50,function(x) mean(num_epiHybWith[x,]/expected[x],na.rm=TRUE))),size=1.5,color=wes_palette("Darjeeling1")[1]) + 
  geom_line(aes(x=2:K,y=sapply(2:50,function(x) mean(num_epiHybBetw[x,]/expected2[x],na.rm=TRUE))),size=1.5,color=wes_palette("Darjeeling1")[2]) +
  geom_hline(aes(yintercept=net_density))


#Proportion of substitutions that can introgress
uppers4 = apply(frac_pos,1,function(x) quantile(x,0.975,na.rm=TRUE))
lowers4 = apply(frac_pos,1,function(x) quantile(x,0.025,na.rm=TRUE))


plot4 = ggplot()+ geom_smooth(data=mpos,aes(x=Var1,y=value))+
  geom_ribbon(aes(x=1:50,ymin=lowers4,ymax=uppers4),fill="grey60",alpha=0.5)+
  coord_cartesian(ylim=c(0,1))+
  labs(x="Substitutions (per population)",y="Fraction of substitutions that can introgress")
