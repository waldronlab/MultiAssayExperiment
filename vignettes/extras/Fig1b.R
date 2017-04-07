library(MultiAssayExperiment)
library(ggplot2)
library(scales)

if(!file.exists("coadMAEO.rds")){
  download.file("http://s3.amazonaws.com/multiassayexperiments/coadMAEO.rds",
                destfile="coadMAEO.rds")
}
if(!file.exists("brcaMAEO.rds")){
  download.file("http://s3.amazonaws.com/multiassayexperiments/brcaMAEO.rds",
                destfile="brcaMAEO.rds")
}

coad <- readRDS("coadMAEO.rds")
brca <- readRDS("brcaMAEO.rds")

fun <- function(mae){
  mae <- mae[, , c("Mutations", "gistict")]
  mae <- mae[, !rownames(colData(mae)) %in% names(duplicated(mae)$Mutations), ]
  mae <- mae[, complete.cases(mae), ]
  mae$nummuts <- sapply(experiments(mae)$Mutations, length) / 50
  mae$cnload <- colMeans(abs(assay(mae[, , "gistict"])[[1]]))
  res <- cbind(mae$nummuts, mae$cnload)
  colnames(res) <- c("N of Mutations", "SCNA level")
  res <- data.frame(res)
  fit <- lm(SCNA.level ~ log10(N.of.Mutations), data=res)
  res$predictions = predict(fit)
  output <- list(res=res, fit=fit)
  return(output)
}

cores <- fun(coad)
brres <- fun(brca)

allRes <- rbind(
    cbind(Cancer="Colon Adenocarcinoma", cores$res),
    cbind(Cancer="Breast", brres$res))
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

png("scna.png", width=4, height=4, units="in", res=600)
ggplot(allRes, aes(N.of.Mutations, SCNA.level, color=Cancer))+
    geom_point()+
    scale_x_continuous(trans = log2_trans(), limits=c(0.2,NA),
        breaks = trans_breaks("log2", function(x) 2^x))+
    geom_smooth(method='lm')+
    theme_classic(15)+
    xlab("Mutations/Mb")+
    ylab("SCNA Score")+
    labs(color="")+
    theme(legend.position="bottom")+
    scale_colour_manual(values=cbPalette)
dev.off()
