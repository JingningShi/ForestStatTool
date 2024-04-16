#' @title 样地分割

#' @author Chaofan Zhou: zhouchaofan@ifrit.ac.cn
#' @author Xiao He: hexiao@ifrit.ac.cn
#' @author Guangshuang Duan: oliverdgs@163.com
#' @author Jingning Shi: shijingning@hebau.edu.cn
#' @author Xiangdong Lei: xdlei@ifrit.ac.cn

#' @usage Plot_Divide(Data, Plot, X, Y, Num_xy = NULL, Length_xy = NULL, Origin = c(0,0), Range_xy)

#' @description 此函数实现了将矩形样地划分为面积相等的若干小样地的功能。可按照坐标X和Y方向上设置的数量或长度划分小样地。

#' @param Data （必须）数据框。需要分割的样地数据框。
#' @param Plot （必须）向量，建议为字符型。样地号。
#' @param X （必须）数值向量。样地林木横坐标。单位：m。
#' @param Y （必须）数值向量。样地林木纵坐标。单位：m。
#' @param Num_xy （可选）数值向量，长度2。均大于1。横纵坐标方向划分小样地的个数。
#' @param Length_xy （可选）数值向量，长度2。均大于0且小于横纵坐标长度的一半。横纵坐标方向划分小样地的长度。单位：m。
#' @param Origin 数值向量，长度2。样地原点的坐标c(X0, Y0)。默认为c(0, 0)。
#' @param Range_xy （必须）数值向量，长度2。均大于0。矩形样地横坐标X和纵坐标Y方向的长度。单位：m。

#' @details 1）可对多个样地组成的数据框进行样地分割的批量处理。
#' @details 2）参数Length_X和Length_Y决定划分的小样地的大小，参数Num_X和Num_Y决定划分的小样地的数量。两组参数可以根据Length_X×Num_X = Xrange和Length_Y×Num_Y = Yrange相互推导，可以两组同时输入，或至少输入一组。当两组同时输入时，Length_X×Num_X和Length_Y×Num_Y分别需要小于或等于Xrange和Yrange，否则，将以参数Num_X和Num_Y为主，推导出Length_X和Length_Y。
#' @details 3）小样地的编号规则：
#' @details 当横坐标和纵坐标方向划分的小样地个数均小于等于10个时，采用列号（0-9）+行号（0-9）的2位数小样地编号。如34表示从左下角开始的第4列第5行的小样地。
#' @details 当横坐标或纵坐标方向划分的小样地个数超过10个时，采用列号（00-99）+行号（00-99）的4位数小样地编号。如0706表示从左下角开始的第8列第7行的小样地。

#' @return 1）输出结果为一个数据框。
#' @return 包括输入的Data数据框的所有列，以及函数产生的新列。
#' @return 2）新输出变量名：
#' @return subplot：分割后的新样地名，格式为“原样地名-小样地编号”。
#' @return subX：小样地的新X坐标，范围（0，Length_X）。单位：m。
#' @return subY：小样地的新X坐标，范围（0，Length_Y）。单位：m。

#' @export Plot_Divide
#' @name Plot_Divide


#' @importFrom deldir duplicatedxy
#' @importFrom dplyr arrange
#' @importFrom dplyr count
#' @importFrom dplyr desc
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr slice_max
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom plyr create_progress_bar
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr "%>%"

#' @examples ## 加载内置数据
#' @examples data(ForestStatTool)
#' @examples #bigplot是一个100m×100m的大样地数据框，用于分割小样地
#' @examples #横纵坐标各划分5份，结果为划分后的25个20m×20m的小样地的合并数据框。
#' @examples #输入方式
#' @examples subplot1 <- Plot_Divide(Data = bigplot, Plot = bigplot$Plot, X = bigplot$X, Y = bigplot$Y,
#' @examples                         Num_xy = c(5,5), Range_xy = c(100,100))
#' @examples head(subplot1)
#' @examples levels(factor(subplot1$subplot))
#' @examples #一样的结果
#' @examples subplot2 <- Plot_Divide(Data = bigplot, Plot = bigplot$Plot, X = bigplot$X, Y = bigplot$Y,
#' @examples                         Length_xy = c(20,20), Range_xy = c(100,100))
#' @examples head(subplot2)
#' @examples levels(factor(subplot2$subplot))
#' @examples #横坐标划分5份纵坐标各划分4份，结果为划分后的16个15m×20m的小样地的合并数据框
#' @examples subplot3 <- Plot_Divide(Data = bigplot, Plot = bigplot$Plot, X = bigplot$X, Y = bigplot$Y,
#' @examples                         Num_xy = c(5,4), Length_xy = c(15,20), Range_xy = c(100,100))
#' @examples head(subplot3)
#' @examples levels(factor(subplot3$subplot))


Plot_Divide <- function(Data, Plot, X, Y, Num_xy = NULL, Length_xy = NULL, Origin = c(0,0), Range_xy){
  
  # 检测dplyr包是否存在，只有不存在时才会安装
  # if(!require ('dplyr')){
    # install.packages("dplyr")
    # #library(dplyr)
  # }
  # 检测tidyverse包是否存在，只有不存在时才会安装
  # if(!require ('tidyverse')){
    # install.packages("tidyverse")
    # #library(tidyverse)
  # }

  if(length(Range_xy)==2){
    RangeX <- Range_xy[1]
    RangeY <- Range_xy[2]
  }else{
    stop("The length of Range_xy should be 2.")
  }
  if(is.null(Num_xy)){
    #未输入Num_xy
    if(is.null(Length_xy)){
      #未输入Length_xy
      stop("Missing parameters of Num_xy or Length_xy")
    }else{
      #输入了Length_xy
      LengthX <- Length_xy[1]
      LengthY <- Length_xy[2]
      if(LengthX > RangeX/2 | LengthX <= 0 | LengthY > RangeY/2 | LengthY <= 0){
        stop("Length should less than one half of Range, respectively.")}
      NumX <- floor(RangeX/LengthX)
      NumY <- floor(RangeY/LengthY)
    }
  }else{
    #输入了Num_xy
    NumX <- Num_xy[1]
    NumY <- Num_xy[2]
    if(NumX < 2 | NumY < 2){
      stop("NumX and NumY should larger than 1.")
    }
    if(is.null(Length_xy)){
      #未输入Length_xy
      LengthX <- RangeX/NumX
      LengthY <- RangeY/NumY
    }else{
      #输入了Length_xy
      LengthX <- Length_xy[1]
      LengthY <- Length_xy[2]
      if(NumX * LengthX > RangeX){
        LengthX <- RangeX/NumX
      }
      if(NumY * LengthY > RangeY){
        LengthY <- RangeY/NumY
      }
    }
  }

  Data$subplot <- Plot
  Data$subX <- X
  Data$subY <- Y
  Data <- na.omit(Data,cols=c("subplot","subX", "subY"))
  plotname<-unique(Data$subplot)
  numx0 <- NumX-1
  numy0 <- NumY-1
  plotall<-data.frame()
  for (i in 1:length(plotname)) {
    datax <- dplyr::filter(Data,subplot==plotname[i])
    for(j in 0:numy0){
      Ylow <- Origin[2]+j*LengthY
      Yhigh <- Origin[2]+(j+1)*LengthY
      for (k in 0:numx0) {
        Xlow <- Origin[1]+k*LengthX
        Xhigh <- Origin[1]+(k+1)*LengthX
        data_new <- dplyr::filter(datax, subY >= Ylow, subY < Yhigh, subX >= Xlow, subX < Xhigh)
        if(nrow(data_new)>0){
          k0 <- k
          j0 <- j
          if(numx0 < 10 & numy0 < 10){
            data_new$subplot <- paste0(data_new$subplot,"-",k0, j0)
          }else{
            if(k<10){k0 <- paste0(0,k)}
            if(j<10){j0 <- paste0(0,j)}
            data_new$subplot <- paste0(data_new$subplot,"-",k0, j0)
          }
          data_new$subX <- data_new$subX-Xlow
          data_new$subY <- data_new$subY-Ylow
          plotall <- rbind(plotall,data_new)
        }
      }
    }
  }
  return(plotall)
}
