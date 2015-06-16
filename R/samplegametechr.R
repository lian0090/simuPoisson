samplegametechr=function(fpchrmat,mapchr,...){
	nmar=nrow(mapchr)
	gamete=vector(length=nmar)
	chrlength=max(mapchr$pos)
	#minimum and maximum number of recombinations cut off based on Morgan map distance
	minrecmb=qpois(0.025,chrlength/100,lower.tail=T)
	maxrecmb=qpois(0.975,chrlength/100,lower.tail=T)
#Parent geno must not be a matrix intead of a data.frame, or it will cause problems.
    gameteOrigin=sample(c(1:2),1)
	gamete=fpchrmat[gameteOrigin,]
	nrecomb=rpois(1,lambda=chrlength/100);while(nrecomb<minrecmb | nrecomb >maxrecmb){nrecomb=rpois(1,lambda=chrlength/100)}
	if(nrecomb>0){
     interval=matrix(0,nrow=1,ncol=2)
     colnames(interval)=c("left","right")
     interval[1,"left"]=0
     interval[1,"right"]=chrlength
     recomb=vector()
  	for(recomb.id in 1:nrecomb){
  		#After each recombination, the new recombination will not occure within 20cM of the previous crossovers. So, the total length of chromomes that a new crossover can occure reduces.
  		#Intead of random sampling cross-over positions along the whole chromosome, the chromosome were divided into different available intervals after each cross-over occured. Each cross break up the original large interval into two smaller intervals, the new smaller intervals are the original large interval that exclude the 20 cM segment where cross-over are not allowed to occur. 
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
changepoints=apply(recomb,1,function(a){which(mapchr$pos>=a)[1]})#the first marker after the recombination event.
newchunkstart=changepoints[seq(1,nrecomb,by=2)]      #the new chunk produced by recombination, only the odd number of recombinations started the new chunk from the other parent, even number of recombinations end the new chunk from the other parent started by the previous odd recombination. For example, a gamete is P1 gamte, 1st recombination along the chromosome starts a new chunk with P2, 2nd recombination end this new chunk of P2, the remaining segments to the end of chromosome is still from P1.
if(nrecomb>1){newchunkend=changepoints[seq(2,nrecomb,by=2)]
if(nrecomb%%2==1){newchunkend=c(newchunkend,nmar)}}else newchunkend=nmar
for(start.id in 1:length(newchunkstart)){
	newchunkinterval=newchunkstart[start.id]:newchunkend[start.id]
	gamete[newchunkinterval]=fpchrmat[-gameteOrigin,newchunkinterval]}
}
return(gamete)		
}