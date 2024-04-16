#' @title 重复坐标点处理1

#' @author Chaofan Zhou: cfzhou2021@163.com
#' @author Xiao He: hexiao@ifrit.ac.cn
#' @author Guangshuang Duan: oliverdgs@163.com
#' @author Jingning Shi: shijingning@hebau.edu.cn
#' @author Xiangdong Lei: xdlei@ifrit.ac.cn

#' @usage Coord_Move(Data, Plot, X, Y, Origin = c(0,0), Range_xy)

#' @description Coord_Move()通过给予重复坐标一个随机的抖动来解决坐标重复问题。

#' @details 1）Coord_Move()会产生新列newX和newY。
#' @details 2）Coord_Move()中得到的新坐标与重复坐标的距离在0.1~0.3m之间，相对于重复坐标的方位是随机的。
#' @details 3）当各样地均不存在重复坐标时，会输出“There are no duplicated coordinates!”的提醒，此时将不会再输出新数据框。
#' @details 4）可对多个样地组成的数据框进行重复标签批量处理。

#' @param Data （必须）数据框。存在坐标重复的样地林木数据框。
#' @param Plot （必须）向量，建议为字符型。样地号。
#' @param X （必须）数字向量。样地林木横坐标。单位：m。
#' @param Y （必须）数字向量。样地林木纵坐标。单位：m。
#' @param Origin 数字向量，长度2。样地原点的坐标c(X0, Y0)。默认为c(0, 0)。单位：m。
#' @param Range_xy （必须）数值向量，长度2，均大于0。矩形样地横坐标X和纵坐标Y方向的长度。单位m。

#' @return 1）输出结果为一个数据框。
#' @return 包括输入的Data数据框的所有列，以及函数产生的新列。
#' @return 2）新输出变量名：
#' @return newX：调整后的新X坐标。
#' @return newY：调整后的新Y坐标。

#' @export Coord_Move
#' @name Coord_Move

#' @import deldir
#' @import stats

#' @examples ## 加载内置数据
#' @examples data(ForestStatTool)
#' @examples ##rawdata是一个存在重复坐标的样地原始数据框
#' @examples ##Coord_Move()通过给予重复坐标一个随机的抖动来解决坐标重复问题。
#' @examples newdata1 <- Coord_Move(Data = rawdata, Plot = rawdata$Plot, X=rawdata$X, Y=rawdata$Y, Range_xy = c(10,10))
#' @examples newdata1
#' @examples ##检查newdata1中的newX和newY中是否还存在重复坐标
#' @examples Coord_Move(Data = newdata1, Plot = newdata1$Plot, X=newdata1$newX, Y=newdata1$newY, Range_xy = c(10,10))
#' @examples #输出“There are no duplicated coordinates!”提醒，说明所有样地中的样木都不存在重复坐标


Coord_Move <- function(Data, Plot, X, Y, Origin = c(0,0),Range_xy){
  # if(!require("deldir")){
    # install.packages("deldir")
    # library(deldir)
  # }
  #主要变量合理性检测
  if(length(Range_xy)==2){
    RangeX <- Range_xy[1]
    RangeY <- Range_xy[2]
  }else{
    stop("The length of Range_xy should be 2.")
  }
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
  #
  Data$newX <- X
  Data$newY <- Y
  plotname<-unique(Plot)
  plotall<-data.frame()
  checkxy <- c()
  for (i in 1:length(plotname)) {
    plotx <- subset(Data,Plot==plotname[i])
    Xi <- X[Plot==plotname[i]]
    Yi <- Y[Plot==plotname[i]]
    #重复坐标点检查
    dupxy <- deldir::duplicatedxy(Xi,Yi)
    # 何潇修改判断写法，2024-4-10
    if( any(dupxy) ){
      checkxy[i] <- T
      for (j in 1:nrow(plotx)) {
        if(dupxy[j]=="TRUE"){
          repeat{
            Angle <- runif(1, min=0, max=pi*2)
            dist <- runif(1, min=0.1, max=0.3)
            plotx$newX[j] <- round(dist*sin(Angle)+Xi[j],2)
            plotx$newY[j] <- round(dist*cos(Angle)+Yi[j],2)
            if(plotx$newX[j]>Origin[1]&plotx$newX[j]<Origin[1]+RangeX&
               plotx$newY[j]>Origin[2]&plotx$newY[j]<Origin[2]+RangeY){
              break
            }
          }
        }
      }
    }else{
      checkxy[i] <- F
    }
    plotall <- rbind(plotall,plotx)
  }
  # 何潇修改判断写法，2024-4-10
  if(!any(checkxy)){
    print("There are no duplicated coordinates!")
  }else{
    print("Duplicated coordinates found! Column of 'newX' and 'newY' were provided in output dataframe.")
    return(plotall)
  }
}
