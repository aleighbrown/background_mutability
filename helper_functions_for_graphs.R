library(ggplot2)
library(cowplot)
library(eulerr)
library(ggpubr)
library(gridExtra)
library(data.table)
library(dgof)
library(grid)
library(pROC)
library(PRROC)
library(gtools)
library(png)
library(ROCR)
library(stringr)
library(dunn.test)
library(ggalluvial)
library(jsonlite)
library(data.table)
library(Biostrings)
library(stringr)
library(stringdist)
library(EnsDb.Hsapiens.v86)
library(rtracklayer)
library(RColorBrewer)

return_single_aa <- function(thre){
  if(any(AMINO_ACID_CODE == thre)){
    return(names(AMINO_ACID_CODE[AMINO_ACID_CODE == thre]))
  }else{
    return("na")
  }
}

make_score_table <- function(dtab, scores = c("fathmmScore","condel","mutability","CountNew","LR","binomial"),negs = c(-1,1,-1,1,1,-1), roundto = 2){
  combinescores <- data.table(scores)
  combinescores[,AUCRoc := rep(0,length(scores))]
  combinescores[,AUCPR := rep(0,length(scores))]
  combinescores[,MCC := rep(0,length(scores))]
  combinescores[,Sens := rep(0,length(scores))]
  setnames(combinescores,"scores","scoreName")
  
  for(i in 1:length(scores)){
    
    combinescores[scoreName == scores[i], MCC := round(get_mat_better(dtab,scores[i],neg = negs[i]), roundto)]
    
    combinescores[scoreName == scores[i], AUCRoc := round(get_auc_better(dtab,scores[i],neg = negs[i]),roundto)]
    
    combinescores[scoreName == scores[i], AUCPR := round(get_pr_better(dtab,scores[i],neg = negs[i]),roundto)]
    
    combinescores[scoreName == scores[i], Sens := round(sens_at10spec(dtab,scores[i]), roundto)]
    
    
  }
  return(combinescores[order(-MCC)])
}

#setwd("/Users/browna6/mahmood/binom")

ci.median2 <- function (x, conf = 0.95) {
  
  n <- nrow(as.matrix(x))
  
  
  L <- qbinom((1 - conf)/2, n, 0.5)
  U <- n - L + 1
  
  order.x <- sort(x)
  res <- list()
  
  if (L >= U | qbinom((1 - conf)/2, n, 0.5) == 0){
    res$ci <- c(median = median(x), lower = median(x), upper = median(x), N = n, min(x),max(x))}
  else{
    res$ci <- c(median = median(x), lower = order.x[L], upper = order.x[n - L + 1], N = n,min(x),max(x))
  }
  
  r <- transpose(data.frame(res$ci))
  names(r) <- c("median","lower","upper","n samples",'minimum',"maximum")

  return(r)
  
}


#takes a numeric decimal and returns a character with the number of decimals you want
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

r2 <- function(x){  
  SSe <- sum((x$resid)^2);  
  observed <- x$resid+x$fitted;  
  SSt <- sum((observed-mean(observed))^2);  
  value <- 1-SSe/SSt;  
  return(value);  
}  

#
assign.subtype <- function(dtab, freqmeasure = "observed_wgs", group_name = "group", table_flag = "n"){
  if(table_flag == "n"){
    if ("outcome" %in% names(dtab)) {
      dtab[outcome == "S", SubType := "Silent"]
      dtab[outcome == "M", SubType := "Missense"]
      dtab[outcome == "N", SubType := "Nonsense"]
      dtab <- dtab[!is.na(SubType)]
    }  
    
    if ("wildtype" %in% names(dtab)){
      dtab[wildtype == mutant & wildtype != "*", SubType := "Silent"]
      dtab[wildtype != mutant & mutant != "*" & wildtype != "*", SubType := "Missense"]
      dtab[mutant == "*" & wildtype != "*", SubType := "Nonsense"]
      dtab <- dtab[!is.na(SubType)]
    } 
  }else{
    dtab[mutAA == AminoAcid & AminoAcid != "*", SubType := "Silent"]
    dtab[AminoAcid != mutAA & mutAA != "*" & AminoAcid != "*", SubType := "Missense"]
    dtab[mutAA == "*" & AminoAcid != "*", SubType := "Nonsense"]
    dtab <- dtab[!is.na(SubType)]
  }




    
    print("HERE, assigned group")

    dtab[is.na(get(freqmeasure)),(group_name) := "0"]
    dtab[get(freqmeasure) == 0,(group_name) := "0"]
    dtab[get(freqmeasure) == 1,(group_name) := "1"]
    dtab[get(freqmeasure) == 2,(group_name) := "2"]
    dtab[get(freqmeasure) >= 3,(group_name) := "3+"]

  
  if("CodingKey" %in% names(dtab)){
    dtab[CodingKey == "Codon", real := "y"]
    dtab[CodingKey != "Codon", real := "n"]
  }
  return(dtab)
}

assign_freqgroup <- function(dtab, freqmeasure = "observed_wgs", group_name = "group"){
  
  dtab[is.na(get(freqmeasure)),(group_name) := "0"]
  dtab[get(freqmeasure) == 0,(group_name) := "0"]
  dtab[get(freqmeasure) == 1,(group_name) := "1"]
  dtab[get(freqmeasure) == 2,(group_name) := "2"]
  dtab[get(freqmeasure) >= 3,(group_name) := "3+"]
  return(dtab)
}

lm_eqn <- function(x, y){
  m <- lm(y ~ x)
  r2 = specify_decimal(summary(m)$adj.r.squared,3)
  as.character(r2)                 
}

g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}



#create table for roc of a wide format data table with different mutabilities
plot.rocs.diffModels <- function(dtab){
  
  #Ratio of pseudocount over mutability
  
  predFreq <- prediction(dtab[,pseudocount],dtab[,label])
  perfFreq <- performance(predFreq, measure = "tpr", x.measure = "fpr")
  
  to_gg <- as.data.table(cbind(unlist(perfFreq@x.values),unlist(perfFreq@y.values)))
  to_gg[,V3 := "Frequency"]
  print(cols)
  for (i in cols){
    
    predRatio <-prediction(dtab[ , pseudocount / get(i)],dtab[,label])
    perfRatio <- performance(predRatio, measure = "tpr", x.measure = "fpr")
    
    to_gg <- rbind(to_gg,cbind(unlist(perfRatio@x.values),unlist(perfRatio@y.values),rep(paste("MR",i),length(unlist(perfRatio@x.values)))),fill = T)
    
    #inverse mutability 
    predMuta <- prediction(dtab[, 1/ get(i)],dtab[,label])
    perfMuta <- performance(predMuta, measure = "tpr", x.measure = "fpr")
    
    to_gg <- rbind(to_gg,cbind(unlist(perfMuta@x.values),unlist(perfMuta@y.values),rep(i,length(unlist(perfMuta@x.values)))),fill = T)
    
  }
  
  
  
  for (i in grep("Binomial", names(dtab), value=TRUE)){
    
    predRatio <-prediction(dtab[ , -get(i)],dtab[,label])
    perfRatio <- performance(predRatio, measure = "tpr", x.measure = "fpr")
    
    to_gg <- rbind(to_gg,cbind(unlist(perfRatio@x.values),unlist(perfRatio@y.values),rep(i,length(unlist(perfRatio@x.values)))),fill = T)
    
  }
  
  
  #condel performance
  predCondel <- prediction(dtab[,condel],dtab[,label])
  perfCondel <- performance(predCondel, measure = "tpr", x.measure = "fpr")
  
  to_gg <- rbind(to_gg,cbind(unlist(perfCondel@x.values),unlist(perfCondel@y.values),rep("Condel",length(unlist(perfCondel@x.values)))),fill = T)
  
  
  
  
  setnames(to_gg, c("V1","V2","V3"), c("False positive rate", "True positive rate","Predictor"))
  
  to_gg[,`False positive rate` := as.numeric(`False positive rate`)]
  to_gg[,`True positive rate` := as.numeric(`True positive rate`)]
  
  return(to_gg)
  
}


###takes the table returned and plots it with a line for the max mCC
plot_mcc_sens_spec <- function(dtab,measure,line = "y",leg = "n"){
  
  dtab_tab <- create_mcc_sens_spec_table(dtab, measure)

  dtab_plot <- ggplot(dtab_tab, aes(x = Cutoff, y =  Value, color = Measure)) + geom_line(size = 1.5) +
    theme(axis.title.y = element_blank(),legend.position = "bottom",legend.direction = "horizontal",legend.justification = "center")
  if(leg == "n"){
    dtab_plot <- dtab_plot + theme(legend.position = "none")
  }
  
  max_cut <- dtab_tab[Value == dtab_tab[Measure == "MCC",max(Value, na.rm = T)],Cutoff]

  if(line == "y"){
    dtab_plot <- dtab_plot + geom_segment(aes(x = max_cut, y = 1, xend = max_cut, yend = 0),color = "black",linetype=4)
  }
  
 return(dtab_plot)
  
}

###anotate figures with the groups

annotate_groups <- function(figure, xvals,yval,fsize,labels){
  for(i in 1:length(xvals)){
    figure <- figure + annotate("text", x = xvals[i], y = yval, label = labels[i], size = fsize)
  }
  return(figure)
}



#extract the x and y ranges of a plot and returns the result as a vector of xmin, xmax, ymin, ymax
extract_ranges <- function(plotName){
  range_vec <- c(layer_scales(plotName)$x$range$range,layer_scales(plotName)$y$range$range)
  return(range_vec)
}

#takes a plot name and a fraction of the plot size as a decimal, uses extract_ranges to get the ranges, then returns the xmin,max for a grob inside that plot, allows an offset from the x axis
get_xranges <- function(plotName, frac, offset = 0){
  #use extract ranges
  range_vec <- extract_ranges(plotName)
  #find x range
  xr <- range_vec[2] - range_vec[1]
  #offset into the left 
  xmax <- (range_vec[2] - offset)
  xmin <- xmax - (xr * frac)
  return(c(xmin,xmax))
}

get_yranges <- function(plotName,frac,offset = 0){
  #use extract ranges
  range_vec <- extract_ranges(plotName)
  #find y range
  yr <- range_vec[4] - range_vec[3]
  #offset down 
  ymax <- (range_vec[4] - offset)
  ymin <- ymax - (yr * frac)
  return(c(ymin,ymax))
}


##takes an ROC curve and draws a retangle on the ROC curve with the given xmin, xmax, ymin, ymax coords, returns a new ggplot objects
draw.highlight_box <- function(plotName, xmin = 0, xmax = 0.1, ymin = 0, ymax = 0.6, al = 0.3){
  #make a data frame that will be given to geom_rec
  rect <- data.frame(xmin=xmin, xmax= xmax, ymin= ymin, ymax= ymax)
  returnPlot <- plotName + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                                     color="grey20", alpha= al, inherit.aes = FALSE)
  return(returnPlot)
}

#takes an ROC curve and returns the ROC zoomed in
make.zoomed.ROC <- function(plot_dtab, xmin = 0, xmax = 0.1, ymin = 0, ymax = 0.65, ab = T,xbreaks = c(0,0.1),ybreaks = c(0,0.2,0.4,0.6),dots = F){


  
  col_l <- sort(unique(plot_dtab[,"colour"]))
  
  

  returnPlot <- ggplot(plot_dtab, aes(x = x, y = y, color = colour
                                      )) + 
    theme(axis.title = element_blank(),legend.position = "none") + 
    coord_cartesian(xlim = c(xmin,xmax),ylim = c(ymin,ymax)) + 
    geom_line(size = 2) +
    scale_y_continuous(breaks = ybreaks) +
    scale_x_continuous(breaks = xbreaks) +
    scale_color_manual(values = col_l) 

  if(ab){
    returnPlot <- returnPlot + geom_abline(size = 1,linetype = 3) 
  }
  if(dots){
    returnPlot <- returnPlot + annotate("point",x = 0.100763, y = 0.6393305, size = 4.5, color = "black") +
      annotate("point",x = 0.100763, y = 0.6393305, size = 4, color = "#E6C618") + 
      annotate("point",x = 0.009997, y = 0.517, size = 4.5, color = "black") +
     annotate("point",x = 0.009997, y = 0.517, size = 4, color = "#BE5A5A") 
      
  }
  

  return(returnPlot)
}


#makes a combined graph, input is either two scatter plots, or one zooomed ROC probably, frac is the fraction of the graph, offse is the y offset between
#howmany lets you say if you have one or two graphs
make_combined <- function(hst,sct1,sct2 = sct1,frac,offse,howMany = 2){

  if(howMany == 2){
    vals <- c(get_xranges(hst,frac),get_yranges(hst,frac))
    hst <- hst + annotation_custom(grob = ggplotGrob(sct1),
                                   xmin = vals[1], xmax = vals[2], 
                                   ymin = vals[3], ymax = vals[4])
    
    vals <- c(get_xranges(hst,frac),get_yranges(hst,frac,offse))
    
    hst <- hst + annotation_custom(grob = ggplotGrob(sct2),
                                   xmin = vals[1], xmax = vals[2], 
                                   ymin = vals[3], ymax = vals[4])
  }else{
    vals <- c(get_xranges(hst,frac),get_yranges(hst,frac,offse))
    hst <- hst + annotation_custom(grob = ggplotGrob(sct1),
                                   xmin = vals[1], xmax = vals[2], 
                                   ymin = vals[3], ymax = vals[4])
  }

  return(hst)
}



#takes a data table,correlation measure, x,y, font size, and returns a nice grob
make_correlation_grob <- function(dtab, corrTy, xval = 0.03, yval = 0.95,RFont){
  #run a  correlation between mutabiltiy and observed mutations
  
  if(corrTy %in% c("s","p")){
    correl <- dtab[,cor.test(AminoMutabilityPan,AminoCountCosmicPan,method = corrTy)]
    #round it to 2 decimals
    rhoR <- round(correl$estimate,2)
    rho_p <- correl$p.value
    
  }else if(corrTy == "l"){  #if we want to show the Adj R Squared (yes, not a correlation)
    muta_lr <- lm(AminoCountCosmicPan ~ AminoMutabilityPan, data = dtab)

    #this is the pvalue
    rho_p <- summary(muta_lr)$coefficients[2,4]

    #this is the adjusted R Squared
    rhoR <- round(summary(muta_lr)$adj.r.squared,2)
    }
  
  
  #rhoR <- paste0(rhoR,",")
  #if the p value is less than 0.001 it's numerically insignifcant to show that
  if(rho_p < 0.001){
      rho_print <- "p << 0.01"
    }else{
      rho_print <- paste("p =", round(rho_p,3))
    }
    #if the corrlation is significant, bold it
  if(rho_p < 0.01){
    #if the correlation is spearman
    if(corrTy == "s"){
        rho2 <- bquote(bold(rho == .(rhoR) ~ .(rho_print)))
        }else if(corrTy == "p"){
        #if it's pearson correlation
        rho2 <- bquote(r == bold(.(rhoR) ~ .(rho_print)))
      }else if(corrTy == "l"){ #it's the linear regression
        rho2 <- bquote(~R^2 == bold(.(rhoR)))
        }
      
    }else{ #if the correlation is not signficant, no bold
      if(corrTy == "s"){
        rho2 <- bquote(rho == .(rhoR) ~ .(rho_print))
      }else if(corrTy == "p"){
        rho2 <- bquote(r == .(rhoR) ~ .(rho_print))
      }else if(corrTy == "l"){
        rho2 <- bquote(R^2 == .(rhoR))
      }

    }
  #make a lovely little grob tree object
  rhogrob <- grobTree(textGrob(rho2, x = xval,  y = yval, hjust=0,vjust = 1,
                                 gp=gpar(fontsize=RFont)))
  return(rhogrob)
  }

make_small_scatter <- function(gene = "none",subtype = "all",fs = 12,Rvj = 0.1, RFont = 12,dtab,fontface = "bold",frequencyType = "observed_wgs",measures){
  #if the substitution type is not all, then take just a copy of that section
  if(subtype != "all"){
    dt <- copy(dtab[SubType == subtype])
  } else{
    dt <- dtab
  }
  #if we want a specfic gene, create the grob for the name and reduce the table to just that gene
  if(gene != "none"){
    #patch
    setnames(dt, "gene", "Gene")
    dt <- dt[Gene == gene]
  }
  
  
  if(subtype == "Missense"){
    col <-  "#00BFC4"
  }else if(subtype == "Nonsense"){
    col <- "#F8766D"
  } else if(subtype == "Silent"){
    col <- "#7CAE00"
  }else{
    col <- "black"
  }
  brksX <- seq(min(dt$AminoMutabilityPan),max(dt$AminoMutabilityPan),length.out = 2)
  xBrLab <- formatC(brksX, format = "e", digits = 1)
  #xBrLab <- formatC(-1 * brksX, digits = 1)
  
  # if(frequencyType == "observed_wgs"){
  #   addL <- "none"
  #   brks <- seq(min(dt$observed_wgs),max(dt$observed_wgs),length.out = 2)
  #   yBrLab <- brks
  #   upper <- max(dt$observed_wgs)
  if(frequencyType == "observed_wgs"){
      addL <- "reg.line"
      brks <- seq(min(dt$AminoCountCosmicPan),max(dt$AminoCountCosmicPan),length.out = 2)
      yBrLab <- brks
      upper <- max(dt$AminoCountCosmicPan)
  }else{
    addL <- "reg.line"
    brks <- seq(0,max(dt$MutationsPer),length.out = 2)
    yBrLab <- c("0",formatC(brks, format = "e", digits = 1)[2])
    upper <- max(dt$MutationsPer)
  }
  
  sct <- ggscatter(dt, x = "AminoMutabilityPan", y = "AminoCountCosmicPan", color = col, add = addL,conf.int = T,size = 3) + 
    theme(axis.title = element_blank()) + 
    theme(axis.text = element_text(size=12)) +
    scale_y_continuous(breaks = brks,label = yBrLab) + 
    scale_x_continuous(breaks = brksX, label = xBrLab) +
    coord_cartesian(ylim = c(0,upper))

  #multiplier for the r vertical adjustment variable 
  rvjMult = 0
  #check that we don't have a gene
  if(gene != "none"){
    
    geneGrob <- grobTree(textGrob(gene, x=0.03,  y=0.95, hjust=0,vjust = 1,
                                  gp=gpar(fontsize=fs, fontface=fontface)))
    rvjMult <- rvjMult + 1
    sct <- sct + annotation_custom(geneGrob)
  }
  #if we have any measures to draw on the plot
  for(ms in measures){

    mGrob <- make_correlation_grob(dt,corrTy = ms, xval = 0.03, yval = 0.95 - (rvjMult * Rvj),RFont = RFont)
    sct <- sct + annotation_custom(mGrob)
    rvjMult <- rvjMult + 1
    
  }
  return(sct)

}

make_inset_ROC <- function(plotName,al = 0.3, xmin = 0.52, xmax = 1.02, ymin = 0, ymax = .5,zoomXmin = 0,
                           zoomXmax = 0.1, zoomYmin = 0, zoomYmax = 0.6,ab = T, xbreaks = c(0,0.1),ybreaks = c(0,0.2,0.4,0.6), dots = F){
  #make the figures with the highlight box  
  wHighLight <- draw.highlight_box(plotName,al = al,xmin = zoomXmin, xmax = zoomXmax, ymin = zoomYmin,ymax = zoomYmax) 
  #get the data
  d <- ggplot_build(wHighLight)$data[[1]]
  #make the zoomed in version of the graph
  zm <- make.zoomed.ROC(d,xmin = zoomXmin, xmax = zoomXmax, ymin = zoomYmin, ymax = zoomYmax + 0.05, ab = ab,xbreaks = xbreaks, ybreaks = ybreaks,dots = dots)
  #make the grob
  g <- ggplotGrob(zm)
  #make them both together
  all <- wHighLight + annotation_custom(grob = g, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  return(all)
}


###sets something from zero to one
range01 <- function(x){(x-min(x))/(max(x)-min(x))}


#takes a data table, a measure, and whether the measure should be neg or pos, returns MCC
get_mat_better <- function(dtab,measure,neg = 1){
  
  predMeas <- prediction(neg * dtab[, get(measure)],dtab[,label])
  perfMeas <- performance(predMeas,measure = "mat")

  mats <- c(max(unlist(perfMeas@y.values),na.rm = T))
  names(mats) <- measure
  
  return(mats)
  
}

get_mat_max_alpha <- function(dtab,measure,neg = 1){
  
  predMeas <- prediction(neg * dtab[, get(measure)],dtab[,label])
  perfMeas <- performance(predMeas,measure = "mat")
  
  mats <- unlist(perfMeas@x.values)[which.max(unlist(perfMeas@y.values))]
  names(mats) <- measure
  
  return(mats)
  
}
#takes a data table, a measure, and whether the measure should be neg or pos, returns the auc for the ROC
get_auc_better <- function(dtab,measure,neg = 1){
  #inverse mutability 
  predMeas <- prediction(neg * dtab[, get(measure)],dtab[,label])
  perfMeas <- performance(predMeas,measure = "auc")
  
  aucs <- c(max(unlist(perfMeas@y.values),na.rm = T))
  names(aucs) <- measure
  
  return(aucs)
  
}

#find the sensitivty at 
sens_at10spec <- function(dtab,measure){
  rocR <- roc(dtab$label, dtab[,get(measure)])
  
  r1 <- coords(rocR,0.9, input = "specificity",ret = c("sensitivity"))
  names(r1) <- measure
  
  return(r1)
}


##using the ‘PRROC’ library to get the auc-pr
get_pr_better <- function(dtab,measure, neg = 1){
  dtab2 <- copy(dtab[!is.na(get(measure))])
  
  meas <- pr.curve(scores.class0 = neg * dtab2[,get(measure)], weights.class0 = dtab2[,label],curve = T)

  aucpr <- meas$auc.integral
  
  names(aucpr) <- measure
  
  return(aucpr)
}



get_prec_recall <- function(dtab,measure,neg = 1){
  #Ratio of pseudocount over mutability
  predMeas <- prediction(neg * dtab[, get(measure)],dtab[,label])
  prefMeas <- performance(predMeas, measure = "prec", x.measure = "rec")
  
  
  to_gg <- as.data.table(cbind(unlist(prefMeas@x.values),
                               unlist(prefMeas@y.values),
                               unlist(prefMeas@alpha.values)))
  to_gg <- cbind(to_gg, rep(measure,nrow(to_gg)))
  
  names(to_gg) <- c("rec", "prec","alpha", "Measure")
  return(to_gg)
}

#create table for roc
get_tpr_fpr <- function(dtab,measure,neg = 1){
  #Ratio of pseudocount over mutability
  predMeas <- prediction(neg * dtab[, get(measure)],dtab[,label])
  prefMeas <- performance(predMeas, measure = "tpr", x.measure = "fpr")

  
  to_gg <- as.data.table(cbind(unlist(prefMeas@x.values),
                               unlist(prefMeas@y.values),
                               unlist(prefMeas@alpha.values)))
  to_gg <- cbind(to_gg, rep(measure,nrow(to_gg)))
  names(to_gg) <- c("fpr", "tpr", "alpha","Measure")
  return(to_gg)
}

#create table for plotting the Specificity, Sensitivity, and Matthews Correlation
get_mcc_sens_spec <- function(dtab,measure,neg = 1){
  #Ratio of pseudocount over mutability
  predMeas <- prediction(neg * dtab[, get(measure)],dtab[,label])
  prefMeas <- performance(predMeas, measure = "sens", x.measure = "spec")
  
  predMeas2 <- prediction(neg * dtab[, get(measure)],dtab[,label])
  prefMeas2 <- performance(predMeas, measure = "mat")
  
  to_gg <- as.data.table(cbind(unlist(prefMeas@x.values),
                               unlist(prefMeas@y.values),
                               unlist(prefMeas@alpha.values),
                               unlist(prefMeas2@y.values)))
  to_gg <- cbind(to_gg, rep(measure,nrow(to_gg)))
  names(to_gg) <- c("spec", "sens", "alpha","MCC","Measure")
  return(to_gg)
}






make_roc_table <- function(dtab, scores = c("fathmmScore","condel","mutability","observed_wgs","MR","binomial"),negs = c(-1,1,-1,1,1,-1)){
  full_table <- data.table()
  for(i in 1:length(scores)){
    to_gg <- get_tpr_fpr(dtab,scores[i],negs[i])
    full_table <- rbind(full_table, to_gg)
  }
  return(full_table)
}

make_pr_table <- function(dtab, scores = c("fathmmScore","condel","mutability","observed_wgs","MR","binomial"),negs = c(-1,1,-1,1,1,-1)){
  full_table <- data.table()
  for(i in 1:length(scores)){
    to_gg <- get_prec_recall(dtab,scores[i],negs[i])
    full_table <- rbind(full_table, to_gg)
  }
  return(full_table)
}

plot_roc_table <- function(full_table){
  plt <- ggplot(full_table, aes(x = fpr, y = tpr, group = Measure)) + geom_line(aes(color = Measure),size = 1.3) + 
    theme(axis.title = element_blank()) +
    scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Dark2"))(9))+
    theme(aspect.ratio=1) 
  return(plt)
}


plot_pr_table <- function(full_table){
  plt <-ggplot(full_table, aes(x = rec, y = prec, group = Measure)) + geom_line(aes(color = Measure),size = 1.3) + 
    theme(axis.title = element_blank()) +
    scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Dark2"))(9))+
    theme(aspect.ratio=1) +  theme(plot.margin=unit(c(0,0,0,1.2),"cm"))
  return(plt)
}

give.n <- function(x){
  return(c(y = median(x) - 0.02, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}


label_vector <- function(smal_table){
  lab_vec <- character()
  for(i in 1:nrow(smal_table)){
    if(smal_table[i,BinnedRole] == "TSG"){
      lab_vec <- c(lab_vec,bquote(bolditalic(.(smal_table[i,Gene]))))
    } else if(smal_table[i,BinnedRole] == "Oncogene"){
      lab_vec <- c(lab_vec,bquote(underline(.(smal_table[i,Gene]))))
    }else{
      lab_vec <- c(lab_vec,smal_table[i,Gene])
    }
  }
  return(lab_vec)
}



make_standardized_values_nuc_mut <- function(gene_nucleotide){
  #get the mutability and counts by nucleotide along the gene
  #gene_nucleotide <- get.geneNucMutability(gene)
  #get the frequency not counting zeros
  stand_obs <- (gene_nucleotide[CosmicPanCount != 0,c("mutabilityPan","context","wildtype", "mutant")])

  stand_obs[,Measure := "Observed"]
  #get the mutability
  stand_all <- (gene_nucleotide[,c("mutabilityPan","context","wildtype","mutant")])
  stand_all[,Measure := "All Mutations"]
  #combine them into a long table form
  standard_table <- rbind(stand_obs,stand_all)
  #flip so the pyramdine
  standard_table[wildtype %in% c("A", "G"),context := as.character(reverseComplement(DNAStringSet(context)))]
  standard_table[wildtype %in% c("A", "G"),mutant := as.character(reverseComplement(DNAStringSet(mutant)))]
  standard_table[wildtype %in% c("A", "G"),wildtype := as.character(reverseComplement(DNAStringSet(wildtype)))]
  standard_table[,PlotMutation := paste0(substr(context,1,1),"[",wildtype,">",mutant,"]",substr(context,3,3))]
  standard_table[, Count := .N, by = c("PlotMutation","Measure")]
  standard_table[,Prop := Count / .N, by = Measure]
  setorder(standard_table,-mutabilityPan)
  setnames(standard_table, "mutabilityPan", "mutability")
  
  return((standard_table))
}

get_negs_vector <-  function(dtab,scores){
  negs = rep(1,length(scores))
  for(s in 1:length(scores)){
    neg = get_auc_better(dtab, scores[s], neg = -1)
    pos = get_auc_better(dtab, scores[s], neg = 1)
    if(neg > pos){
      negs[s] = -1
    }else{
      negs[s] = 1
    }
  }
  return(negs)
}
