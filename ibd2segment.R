#tract.file should be beagle fibd or ibd tracts file. This link might help if you don't know how to get these files (http://faculty.washington.edu/sguy/beagle/ibd_and_hbd/ibd_and_hbd.html).
#marker.file should be in beagle markers file format (rsID Centimorgan Allele_1 Allele_2, only first two columns are used).
#bim.file should be in plink bim file format (actually only first 4 columns are needed, "Chr rsID Centimorgan base_pair_position")
#fam.file should be in plink fam file format (six columns: FID IID PID MID Gender Affection_status), Affection_status can't be all missing if you want to run plink shared segment analysis later.

ibd2segment <- function(tract.file, marker.file, bim.file, fam.file, out.file, segment.file, fibd=TRUE){
	ibd <- read.table(ibd.file,as.is=TRUE)[,1:4]
	names(ibd) <- c("IID1","IID2","StartIndex","EndIndex")
	fam <- read.table(fam.file,as.is=TRUE)[,c(1,2,6)]
	names(fam) <- c("FID","IID","Aff")
	tem <- merge(ibd,fam,by.x="IID1",by.y="IID",sort=F,all.x=T)
	names(tem)[grep("^FID$",names(tem))] <- "FID1"
	names(tem)[grep("^Aff$",names(tem))] <- "Aff1"
	tem <- merge(tem,fam,by.x="IID2",by.y="IID",sort=F,all.x=T)
	names(tem)[grep("^FID$",names(tem))] <- "FID2"
	names(tem)[grep("^Aff$",names(tem))] <- "Aff2"
	tem$PHE <- -9
	tem$PHE[tem$Aff1==2 & tem$Aff2==2] <- 1
	tem$PHE[tem$Aff1==1 & tem$Aff2==1] <- -1
	tem$PHE[tem$Aff1==1 & tem$Aff2==2] <- 0
	tem$PHE[tem$Aff1==2 & tem$Aff2==1] <- 0
	marker <- read.table(marker.file,as.is=TRUE)[,1:2]
	names(marker) <- c("SNP","cM")
	bim <- read.table(bim.file,as.is=TRUE)[,c(1,2,4)]
	names(bim) <- c("CHR","SNP","BP")
	snp <- merge(bim,marker,sort=FALSE)
	tem$BP1 <- snp$BP[(tem$StartIndex+1*fibd)]
	tem$BP2 <- snp$BP[(tem$EndIndex+1*fibd)]
	tem$SNP1 <- snp$SNP[(tem$StartIndex+1*fibd)]
	tem$SNP2 <- snp$SNP[(tem$EndIndex+1*fibd)]
	tem$cM1 <- snp$cM[(tem$StartIndex+1*fibd)]
	tem$cM2 <- snp$cM[(tem$EndIndex+1*fibd)]
	tem$NSNP <- tem$EndIndex - tem$StartIndex + 1
	tem$KB <- (tem$BP2 - tem$BP1)/1000
	tem$cM <- tem$cM2 -tem$cM1
	tem$CHR <- snp$CHR[1]
	write.table(tem[,c("FID1","IID1","FID2","IID2","PHE","CHR","BP1","BP2","SNP1","SNP2","NSNP","KB")],file=segment.file,row.names=FALSE,quote=FALSE,sep="\t")
}
