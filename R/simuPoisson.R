simuPoisson<-function(fpgenomat,map,gtype,nptotal,...){
savedir=list(...)$savedir
project=list(...)$project
fpgenomat=as.matrix(fpgenomat)
if(nrow(map)!=ncol(fpgenomat)){stop("fpgenomat must has the same number of markers as map")}
#there are a lot of markers whose map positions are exactly the same. 
nmar=nrow(map)
simugeno=matrix(NA,nrow=nptotal,ncol=nmar)
#only simulate poisson distribution for F2 populations.
if(gtype=="F2"){
#set.seed(1)
for(i in 1:nptotal){
gamete1=samplegamete(fpgenomat,map)
gamete2=samplegamete(fpgenomat,map)
simugeno[i,]=(gamete1+gamete2)/2
}
cmbsimugeno=rbind(map$chr,map$pos,fpgenomat,simugeno)
names(cmbsimugeno)=map$mkname
cmbsimugeno=cbind(LINE=c(rownames(cmbsimugeno)[1:4],c(1:nptotal)),cmbsimugeno)
if(!is.null(savedir)){
dir.create(savedir,recursive=T)
saveRDS(cmbsimugeno,file=file.path(savedir,paste(project,".rds",sep="")))	
}
return(cmbsimugeno)
}
}  

			
