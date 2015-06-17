simuPoisson<-function(parentsGeno,chr,cM,N){
  #parentsGeno
  # 2 x p matrix for the parents genotypes. Each marker is a column.
  #chr: a vector of size p for the chromosome numbers of all markers
  #cM: a vector of size p for the centiMorgan map positions for each marker.
  #N: size of the population to simulate
  parentsGeno=as.matrix(parentsGeno)
  if(any(apply(parentsGeno,2,is.character))){
  	stop("parents genotypes must be coded additively: coding for heterozygote is half of two homozygotes genotype ")
  }
  nmar=length(chr)

  if(nmar!=ncol(parentsGeno)){stop("parentsGeno must have the same number of markers as in chr")}
  if(nmar!=length(cM)){stop("chr must be the same length as cM")}

  simugeno=matrix(NA,nrow=N,ncol=nmar)
  
  for(i in 1:N){
    gamete1=samplegamete(parentsGeno,chr,cM)
    gamete2=samplegamete(parentsGeno,chr,cM)
    simugeno[i,]=(gamete1+gamete2)
  }
  return(simugeno)
}


##A function to sample gamete for poisson distribution
samplegamete=function(parentsGeno,chr,cM){
  parentsGeno=as.matrix(parentsGeno)
  chrlevels=unique(map$chr)
  nchr=length(chrlevels)
  gamete=vector(length=length(chr))
  for(i in 1:nchr){
  	chr.id=chrlevels[i]
    index.chr=which(map$chr==chr.id)
    Geno.i=parentsGeno[,index.chr]
    cM.i=cM[index.chr]
    gamete[index.chr]=samplegametechr(Geno.i,cM.i)
  }
  #check whether you got right
  #cat(range(newchunkinterval),"\n",all(gamete[newchunkinterval]==parentsGeno[-gameteOrigin,newchunkinterval],na.rm=T),"\n")}
  return(gamete)
}

##sample gamate for each chromosome

samplegametechr=function(parentsGeno,cM){
  nmar=length(cM)
  gamete=vector(length=nmar)
  chrlength=max(cM)
  #minimum number of recombinations based on Morgan map distance
  minrecmb=qpois(0.025,chrlength/100,lower.tail=T)
  #maximum number of recombinations based on Morgan map distance
  maxrecmb=qpois(0.975,chrlength/100,lower.tail=T)
  gameteOrigin=sample(c(1:2),1)
  gamete=parentsGeno[gameteOrigin,]
  nrecomb=rpois(1,lambda=chrlength/100)
  while(nrecomb<minrecmb | nrecomb >maxrecmb){nrecomb=rpois(1,lambda=chrlength/100)}
  if(nrecomb>0){
    interval=matrix(0,nrow=1,ncol=2)
    colnames(interval)=c("left","right")
    interval[1,"left"]=0
    interval[1,"right"]=chrlength
    recomb=vector()
    for(recomb.id in 1:nrecomb){
      #After each recombination, the new recombination will not occure within 20cM of the previous crossovers. So, the total length of chromomes that a new crossover can occure reduces.
      #Intead of random sampling cross-over positions along the whole chromosome, the chromosome were divided into different available intervals after each cross-over occured. Each cross-over break up the original large interval into two smaller intervals, the new smaller intervals are the original large interval that exclude the 20 cM segment where cross-over are not allowed to occur. 
      ninterval=nrow(interval)
      length.interval=apply(interval,1,function(a){a[2]-a[1]})
      total.length=sum(length.interval)
      occured.interval=sample(1:ninterval,1,prob=length.interval/total.length)
      recomb.i=runif(1,interval[occured.interval,"left"],interval[occured.interval,"right"])
      #new intervals is always the last interval no matter where the breaks occurs. So, the interval orders are not related to their map positions
      #the new break point will not affect the left and right bounder of the occured original larger interval.
      newinterval=ori.large.interval=interval[occured.interval,]
      #the right bound of the new interval inherit directly from the original larger interval
      newinterval["left"]=min(recomb.i+10,ori.large.interval["right"])
      interval[occured.interval,"right"]=max(ori.large.interval['left'],recomb.i-10) #the right bounder might be the same as the left bound, so this interval length will turn to 0 in that case.
      interval=rbind(interval,newinterval)
      #now sort interval according to their map position
      interval=interval[order(interval[,1]),]
      recomb=c(recomb,recomb.i)
    }
    recomb=matrix(sort(recomb))
    changepoints=apply(recomb,1,function(a){which(cM>=a)[1]})#the first marker after the recombination event.
    newchunkstart=changepoints[seq(1,nrecomb,by=2)]      #the new chunk produced by recombination, only the odd number of recombinations started the new chunk from the other parent, even number of recombinations end the new chunk from the other parent started by the previous odd recombination. For example, a gamete is P1 gamte, 1st recombination along the chromosome starts a new chunk with P2, 2nd recombination end this new chunk of P2, the remaining segments to the end of chromosome is still from P1.
    if(nrecomb>1){
    	newchunkend=changepoints[seq(2,nrecomb,by=2)]
       if(nrecomb%%2==1){
       	newchunkend=c(newchunkend,nmar)
       	}
      }else { #only one recombinant
             newchunkend=nmar
        }
      n.newchunk=length(newchunkstart) 
    for(i in 1:n.newchunk){
      newchunkinterval=newchunkstart[i]:newchunkend[i]
      gamete[newchunkinterval]=parentsGeno[-gameteOrigin,newchunkinterval]}
  }
    #gamete genotype is coded as half the value of the inbred parents
  return(gamete/2)		
}



