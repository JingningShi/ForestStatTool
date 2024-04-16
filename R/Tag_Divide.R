#' @title 重复标签处理1

#' @author Chaofan Zhou: cfzhou2021@163.com
#' @author Xiao He: hexiao@ifrit.ac.cn
#' @author Guangshuang Duan: oliverdgs@163.com
#' @author Jingning Shi: shijingning@hebau.edu.cn
#' @author Xiangdong Lei: xdlei@ifrit.ac.cn

#' @usage Tag_Divide(Data, Plot, Tag)

#' @description Tag_Divide通过在重复标签后添加新标注来解决标签重复问题。

#' @details 1）只有Tag_Divide会产生新列newTag。
#' @details 2）当各样地均不存在重复标签时，会输出“There are no duplicated tags!”的提醒，此时将不会再输出新数据框。
#' @details 3）可对多个样地组成的数据框进行重复标签批量处理。

#' @param Data （必须）数据框。存在标签重复的样地林木数据框。
#' @param Plot （必须）向量，建议为字符型。样地号。
#' @param Tag （必须）向量，建议为字符型。样地林木标签。

#' @return 1）输出结果为一个数据框。
#' @return 包括输入的Data数据框的所有列，以及函数可能产生的新列。
#' @return 2）新输出变量名：
#' @return newTag：重新编号的新标签，格式为“原标签-编号”

#' @export Tag_Divide
#' @name Tag_Divide

#' @examples ## 加载内置数据
#' @examples data(ForestStatTool)
#' @examples ##rawdata是一个存在重复标签的样地原始数据框
#' @examples rawdata
#' @examples ##Tag_Divide()通过在重复标签后添加新标注来解决标签重复问题。
#' @examples newdata1 <- Tag_Divide(Data = rawdata, Plot = rawdata$Plot, Tag = rawdata$Tag)
#' @examples newdata1
#' @examples ##检查新列newTag是否还存在重复标签
#' @examples Tag_Divide(Data = newdata1, Plot = newdata1$Plot, Tag = newdata1$newTag)
#' @examples ##输出“There are no duplicated tags!”提醒，说明所有样地中的样木都不存在重复标签


Tag_Divide <- function(Data, Plot, Tag){
  #主要变量合理性检测
  if(sum(is.na(Plot))!=0){
    stop("Missing value (or NA) in 'Plot'.")
  }
  if(sum(is.na(Tag))!=0){
    stop("Missing value (or NA) in 'Tag'.")
  }
  Data$newTag <- Tag
  plotname<-unique(Plot)
  plotall<-data.frame()
  checktags <- c()
  for (i in 1:length(plotname)) {
    plotx <- subset(Data,Plot==plotname[i])
    tags <- plotx$newTag
    #重复标签检查
    duptags <- duplicated(tags)
    # if(T %in% duptags){ }# 何潇修改判断语句写法,2024-4-10
    if( any(duptags) ){
      checktags[i] <- T
      for (j in 1:nrow(plotx)) {
        if(duptags[j]=="TRUE"){
          k <- 1
          repeat{
            plotx$newTag[j] <- c(paste(tags[j],k,sep = "_"))
            if(sum(plotx$newTag[j]==tags)==0){
              break
            }
            k <- k+1
          }
        }
      }
    }else{
      checktags[i] <- F
    }
    plotall <- rbind(plotall,plotx)
  }
  # if(!(T %in% checktags)){ }# 何潇修改判断语句写法,2024-4-10
  if( !any(checktags) ){
    print("There are no duplicated Tags!")
  }else{
    print("Duplicated Tags found! Column of 'newTag' was provided in output dataframe.")
    return(plotall)
  }
}

