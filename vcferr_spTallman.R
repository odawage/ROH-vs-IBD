
# https://github.com/spTallman/vcfErr.git

library(data.table)
library(stringr)
# 
# vcf = "/net/fs-2/scale/OrionStore/Scratch/odwa/paper-1/NE200_100cm_rep0/Half_Dens/inped_2.vcf"
# errProb = 0.01
# outpath = "/net/fs-2/scale/OrionStore/Scratch/odwa/paper-1/NE200_100cm_rep0/Half_Dens"


vcferr <- function(vcf, errProb, outpath) {
   
    vcfFile <- fread(vcf, header=T)
    cat(paste(vcf, "with error rate:", errProb, sep=" "), "\n")
     
    homhetSwitcher <- function(x) {
        
        strand <- sample(c(2,4), replace=F, size=1)
        splitGT <- str_split(x, "|")[[1]]
        
        if (splitGT[strand] == "1") {
            splitGT[strand] <- "0"
            x <- paste(splitGT, collapse="")
        } else if (splitGT[strand] == "0") {
            splitGT[strand] <- "1"
            x <- paste(splitGT, collapse="")
        }
    }
    

    counter=0
    err.counter = 0
    for (i in 10:ncol(vcfFile)) {
        counter <- counter + 1 
        #cat(paste("sample:", colnames(vcfFile)[i], sep=" "), "\n")
        indexErr <- which(rbinom(nrow(vcfFile), length(nrow(vcfFile)), prob=errProb) == 1)
        err.counter <- err.counter + length(indexErr)
        #cat(paste("Number of simulated errors in VCF: ", as.character(length(indexErr))), "\n")
        vcfFile[[i]][indexErr] <- unlist(lapply(as.list(vcfFile[[i]][indexErr]), FUN=function(x2) homhetSwitcher(x2)))
    }
    cat(paste0("The average amount of error generated: ",err.counter/ncol(vcfFile), "\n"))
    vcfPrefix <- as.character(strsplit(vcf, ".vcf")[[1]])
    cat(paste("Writing output to ", vcfPrefix, ".err.vcf", sep=""), "\n")
    write.table(vcfFile, paste0(outpath,errProb,".err.vcf", sep=""), sep="\t", col.names=F, row.names=F, quote=F)

}

# errProb = 0.5
# vcf = "vcfErr/example.vcf"
# vcf = "/mnt/users/odwa/zooroh/Code/Dataset/Bulls_Holstein/bulls_200.vcf"
# vcfFile <- fread(vcf, header=T)
# write.table(vcfFile[1:20,1:13],"mini_np.vcf", sep = "\t",col.names=T, row.names=F, quote = F)


