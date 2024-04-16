#' @title 重复坐标点处理2

#' @author Chaofan Zhou: cfzhou2021@163.com
#' @author Xiao He: hexiao@ifrit.ac.cn
#' @author Guangshuang Duan: oliverdgs@163.com
#' @author Jingning Shi: shijingning@hebau.edu.cn
#' @author Xiangdong Lei: xdlei@ifrit.ac.cn

#' @usage Coord_Remove(Data, Plot, X, Y, D)

#' @description Coord_Remove()通过保留重复坐标中胸径最大的林木来解决坐标重复问题。

#' @details 1）当各样地均不存在重复坐标时，会输出“There are no duplicated coordinates!”的提醒，此时将不会再输出新数据框。
#' @details 2）可对多个样地组成的数据框进行重复标签批量处理。

#' @param Data （必须）数据框。存在坐标重复的样地林木数据框。
#' @param Plot （必须）向量，建议为字符型。样地号。
#' @param X （必须）数字向量。样地林木横坐标。单位：m。
#' @param Y （必须）数字向量。样地林木纵坐标。单位：m。
#' @param D （必须）数字向量。样地林木胸径。单位：cm。

#' @return 1）输出结果为一个数据框。
#' @return 包括输入的Data数据框的所有列，以及函数产生的新列。
#' @return 2）新输出变量名：
#' @return newX：调整后的新X坐标。
#' @return newY：调整后的新Y坐标。

#' @export Coord_Remove
#' @name Coord_Remove

#' @import deldir
#' @import dplyr

#' @examples ## 加载内置数据
#' @examples data(ForestStatTool)
#' @examples ##rawdata是一个存在重复坐标的样地原始数据框
#' @examples ##Coord_Remove()通过保留重复坐标中胸径最大的林木来解决坐标重复问题。
#' @examples ##参数D有空值，删除存在空值的行
#' @examples rawdata0 <- subset(rawdata,!is.na(rawdata$D))
#' @examples newdata2 <- Coord_Remove(Data = rawdata0, Plot = rawdata0$Plot, X = rawdata0$X, Y = rawdata0$Y, D = rawdata0$D)
#' @examples newdata2
#' @examples ##检查newdata2中的newX和newY中是否还存在重复坐标
#' @examples Coord_Remove(Data = newdata2, Plot = newdata2$Plot, X = newdata2$X, Y = newdata2$Y, D = newdata2$D)
#' @examples #输出“There are no duplicated coordinates!”提醒，说明所有样地中的样木都不存在重复坐标


Coord_Remove <- function(Data, Plot, X, Y, D){
  # if(!require("deldir")){
    # install.packages("deldir")
    # library(deldir)
  # }
  # if(!require("dplyr")){
    # install.packages("dplyr")
    # library(dplyr)
  # }
  #主要变量合理性检测
  if(sum(is.na(Plot))!=0){
    stop("Missing value (or NA) in 'Plot'.")
  }
  if(sum(is.na(X))!=0){
    stop("Missing value (or NA) in 'X'.")
  }else if(!is.numeric(X)){
    stop("'X' must be numeric.")
  }
  if(sum(is.na(Y))!=0){
    stop("Missing value (or NA) in 'Y'.")
  }else if(!is.numeric(Y)){
    stop("'Y' must be numeric.")
  }
  if(sum(is.na(D))!=0){
    stop("Missing value (or NA) in 'D'.")
  }else if(!is.numeric(D)){
    stop("'D' must be numeric.")
  }
  #
  Data$newX <- X
  Data$newY <- Y
  Data$DD <- D
  plotname<-unique(Plot)
  plotall<-data.frame()
  checkxy <- c()
  for (i in 1:length(plotname)) {
    plotx <- subset(Data,Plot==plotname[i])
    plotx <- plotx[order(plotx$DD,decreasing = T),]
    dupxy <- deldir::duplicatedxy(plotx$newX,plotx$newY)
    #重复坐标点检查
    # 何潇修改判断条件写法，2024-4-10
    if( any(dupxy) ){
      checkxy[i] <- T
      plotx$dupxy <- dupxy
      plotx <- subset(plotx,dupxy=="FALSE")
      plotall <- rbind(plotall,plotx)
    }else{
      checkxy[i] <- F
      plotx$dupxy <- dupxy
      plotall <- rbind(plotall,plotx)
    }
  }
  # 何潇修改判断条件写法，2024-4-10
  if( !any(checkxy) ){
    print("There are no duplicated coordinates!")
  }else{
    print("Duplicated coordinates found and have been removed.")
    plotall <- plotall[,!colnames(plotall)%in%c("newX","newY","DD","dupxy")]
    plotall <- dplyr::arrange(plotall,Tag)
    return(plotall)
  }
}
