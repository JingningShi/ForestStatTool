#' @title Deal with duplicate tags by DBH. 基于胸径的重复标签处理

#' @author Chaofan Zhou: cfzhou2021@163.com
#' @author Xiao He: hexiao@ifrit.ac.cn
#' @author Guangshuang Duan: oliverdgs@163.com
#' @author Jingning Shi: shijingning@hebau.edu.cn
#' @author Xiangdong Lei: xdlei@ifrit.ac.cn

#' @usage Tag_Remove(Data, Plot, Tag, D)

#' @description Tag_Remove通过删除标签重复的较小胸径的林木而保留标签重复的最大胸径的林木来解决标签重复问题。

#' @details 1）当各样地均不存在重复标签时，会输出“There are no duplicated tags!”的提醒，此时将不会再输出新数据框。
#' @details 2）可对多个样地组成的数据框进行重复标签批量处理。

#' @param Data （必须）数据框。存在标签重复的样地林木数据框。
#' @param Plot （必须）向量，建议为字符型。样地号。
#' @param Tag （必须）向量，建议为字符型。样地林木标签。
#' @param D （必须）数字向量。样地林木胸径。单位：cm。

#' @return 1）输出结果为一个数据框。
#' @return 包括输入的Data数据框的所有列，以及函数可能产生的新列。
#' @return 2）新输出变量名：
#' @return newTag：重新编号的新标签，格式为“原标签-编号”

#' @export Tag_Remove
#' @name Tag_Remove

#' @importFrom dplyr arrange

#' @examples ## 加载内置数据
#' @examples data(ForestStatTool)
#' @examples ##rawdata是一个存在重复标签的样地原始数据框
#' @examples rawdata
#' @examples ##Tag_Remove()通过保留重复标签中胸径最大的林木来解决标签重复问题。
#' @examples ##参数D有空值，删除存在空值的行
#' @examples rawdata0 <- subset(rawdata,!is.na(rawdata$D))
#' @examples newdata2 <- Tag_Remove(Data = rawdata0, Plot = rawdata0$Plot, Tag = rawdata0$Tag, D = rawdata0$D)
#' @examples newdata2
#' @examples ##检查newdata2中的Tag是否还存在重复标签
#' @examples Tag_Remove(Data = newdata2, Plot = newdata2$Plot, Tag = newdata2$Tag, D = newdata2$D)
#' @examples #输出“There are no duplicated tags!”提醒，说明所有样地中的样木都不存在重复标签


Tag_Remove <- function(Data, Plot, Tag, D){
  if(!require("dplyr")){
    install.packages("dplyr")
    library(dplyr)
  }

  #主要变量合理性检测
  if(sum(is.na(Plot))!=0){
    stop("Missing value (or NA) in 'Plot'.")
  }
  if(sum(is.na(Tag))!=0){
    stop("Missing value (or NA) in 'Tag'.")
  }
  if(sum(is.na(D))!=0){
    stop("Missing value (or NA) in 'D'.")
  }else if(!is.numeric(D)){
    stop("'D' must be numeric.")
  }
  Data$newTag <- Tag
  Data$DD <- D
  plotname<-unique(Plot)
  plotall<-data.frame()
  checktags <- c()
  for (i in 1:length(plotname)) {
    plotx <- subset(Data,Plot==plotname[i])
    plotx <- plotx[order(plotx$DD,decreasing = T),]
    tags <- plotx$newTag
    #重复标签检查
    duptags <- duplicated(tags)
    # if(T %in% duptags){ }# 何潇修改判断语句写法,2024-4-10
    if( any(duptags) ){
      checktags[i] <- T
      plotx$duptags <- duptags
      plotx <- subset(plotx,duptags=="FALSE")
    }else{
      checktags[i] <- F
      plotx$duptags <- duptags
    }
    plotall <- rbind(plotall,plotx)
  }
  if(!(T %in% checktags)){
    cat("There are no duplicated Tags!\n")
  }else{
    cat("Duplicated Tags found and have been removed.\n")
    plotall <- plotall[,!colnames(plotall) %in% c("newTag","DD","duptags")]
    plotall <- arrange(plotall,Tag)
    return(plotall)
  }
}




