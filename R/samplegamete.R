##A function to sample gamete for poisson distribution
samplegamete=function(fpgenomat,map=map,...){
fpgenomat=as.matrix(fpgenomat)
map=data.frame(map);map$mkname=as.character(map$mkname)
nchr=length(unique(map$chr))
gamete=vector(length=nrow(map))
for(chr.id in 1:nchr){
	mapchr=subset(map,map$chr==chr.id)
    wh.chr=which(map$chr==chr.id)
    fpchrmat=fpgenomat[,wh.chr]
	gamete[wh.chr]=samplegametechr(fpchrmat,mapchr)
	}
#check whether you got right
#cat(range(newchunkinterval),"\n",all(gamete[newchunkinterval]==fpgenomat[-gameteOrigin,newchunkinterval],na.rm=T),"\n")}
return(gamete)
}
