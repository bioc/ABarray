#- $Id: doVennDiagram.R,v 1.1.1.1 2006/06/06 22:06:37 sunya Exp $

#- Draw Venn diagram

"doVennDiagram" <-
function (a, b, c = NULL, names, ...)
{
  if (  is.null(c) ) {
    list.all <- union(a, b)
    list.mat <- matrix(0, nrow = length(list.all), ncol = 2)
    colnames(list.mat) <- c("list1", "list2")
    for(i in 1:length(list.all)) {
      list.mat[i,1] <- list.all[i] %in% a
      list.mat[i,2] <- list.all[i] %in% b
    }
    list.venn <- vennCounts(list.mat)
    drawVennDiagram(list.venn, names = names, mar = rep(1, 4), cex = 1.5, ...)
    ab <- intersect(which(list.mat[,1] == 1), which(list.mat[,2] == 1))
    list.ab <- vector(length = length(list.all))
    list.ab[ab] <- TRUE
    fileName <- "Venn_list_1+2.csv"
    if(! missing(names)) {
      fileName <- paste("Venn", names[1], names[2], ".csv", sep = "_")
    }
    write.table(list.all[list.ab], file = fileName, sep = ",")
    print(paste("List information is written to file", fileName))
    invisible(list.all[list.ab])
  }
  else {
    list.all <- union(a, union(b, c))
    list.mat <- matrix(0, nrow = length(list.all), ncol = 3)
    colnames(list.mat) <- c("list1", "list2", "list3")
    for(i in 1:length(list.all)) {
      list.mat[i,1] <- list.all[i] %in% a
      list.mat[i,2] <- list.all[i] %in% b
      list.mat[i,3] <- list.all[i] %in% c
    }
    list.venn <- vennCounts(list.mat)
    drawVennDiagram(list.venn, names = names, mar = rep(1, 4), cex = 1.5, ...)
    ab <- intersect(which(list.mat[,1] == 1), which(list.mat[,2] == 1))
    ac <- intersect(which(list.mat[,1] == 1), which(list.mat[,3] == 1))
    bc <- intersect(which(list.mat[,2] == 1), which(list.mat[,3] == 1))
    abc <- intersect(which(list.mat[,1] == 1), bc)
    vd <- matrix(0, nrow = length(list.all), ncol = 4)
    rownames(vd) <- list.all
    vd[,1][ab] <- 1
    vd[,2][ac] <- 1
    vd[,3][bc] <- 1
    vd[,4][abc] <- 1
    colnames(vd) <- c("list 1_2", "list 1_3", "list 2_3", "list 1_2_3")
    if(!missing(names)) {
      colnames(vd) <- c(paste(names[1], names[2], sep = "+"),
                        paste(names[1], names[3], sep = "+"),
                        paste(names[2], names[3], sep = "+"),
                        paste(names[1], names[2], names[3], sep = "+"))
    }
    fileName <- "Venn_list_1+2+3.csv"
    if(! missing(names)) {
      fileName <- paste("Venn", names[1], names[2], names[3], ".csv", sep = "_")
    }
    write.table(vd, file = fileName, sep = ",", col.names = NA)
    print(paste("List information is written to file", fileName))
    invisible(vd)
  }
}

######################################################
#- $Log: doVennDiagram.R,v $
#- Revision 1.1.1.1  2006/06/06 22:06:37  sunya
#- ABarray project converted from ab1700 project
#-
#- Revision 1.3  2006/03/14 19:48:30  sunya
#- Changed icp (internal control probe) QC plots.
#- Added function for icp -> icpPlot
#- ANOVA analysis now performs probe filtering, but no FDR is calculated.
#- hclusterPlot now calculate correlation coefficient for probes, previously
#- it used Euclidian distance. The distance between arrays is still Euclidean.
#-
