#' @title 与距离有关的竞争指数_基础版_3个

#' @author Guangshuang Duan: oliverdgs@163.com
#' @author Chaofan Zhou: cfzhou2021@163.com
#' @author Xiao He: hexiao@ifrit.ac.cn
#' @author Mengli Zhou: 1257304343@qq.com
#' @author Jingning Shi: shijingning@hebau.edu.cn
#' @author Xiangdong Lei: xdlei@ifrit.ac.cn

#' @usage CI.dist.base(Data = NULL, Plot, Tag, X, Y, D, SP = NULL, H = NULL, CR = NULL, Neighbor = "number", k = 4, Search_dist = NULL, Shape = "rectangle", Correct = "single", Origin = c(0,0), Range_xy = NULL, Radius = NULL, Buf_dist = 5, Bind = FALSE)

#' @description 此函数实现了利用单木胸径和位置计算与距离有关的单木竞争指数，包括与距离有关的3个变量，每个指标均可分为总竞争、种内竞争和种间竞争。
#' @description 函数对矩形样地和圆形样地均适用（参数Shape），并内置3种边缘矫正方法（参数Correct）。

#' @details 1）输入数据可以选择两种方法：第一种是按照参数Plot, Tag, X, Y, D, SP逐个输入其中SP可选择不输入；第二种是首先修改数据的列名，包括Plot, Tag, X, Y, D, SP，其中SP可不包含在数据中，然后通过Data参数一次性输入。
#' @details 2）当参数Shape="rectangle"（即方形样地）时，3种边缘矫正方法均可用，当参数Shape ="circle "时（即圆形样地时），Correct只能选择逐株优化法（"single"）或固定距离的缓冲区法（"buffer"）。当参数Neighbor="number"时（即固定数量法）时，3种边缘矫正方法均可用，当参数Neighbor="distance"时（即固定半径法）时，边缘矫正方法只能选择为平移法（"translation"）和固定距离的缓冲区法（"buffer"），这时缓冲区距离Buf_dist就等于固定半径法的半径Search_dist。
#' @details 3）如果SP不输入则只计算总竞争，输入SP后会同时计算3个竞争指数的总竞争和种内种间竞争。
#' @details 4）函数内置运算前检查功能，包括判断各个输入参数是否满足计算要求，是否存在重复坐标和标签。

#' @param Data （可选）数据框。样地的所有林木属性信息，当Data输入时，Data中必须包含名称为Plot, Tag, X, Y, D。在Bind = TRUE 时，Data必须输入。
#' @param Plot （必须）向量，建议为字符型。样地号。可以是数值型，也可以是字符型。
#' @param Tag （必须）向量。林木标签（编号）。可以是数值型，也可以是字符型。
#' @param X （必须）数字向量。样地林木横坐标。单位：m。
#' @param Y （必须）数字向量。样地林木纵坐标。单位：m。
#' @param D （必须）数字向量。样地林木胸径。单位：cm。
#' @param SP （可选）向量，可以是数字型，也可以是字符型。树种。
#' @param H （可选）数值向量。样地林木树高。单位：m。
#' @param CR （可选）数值向量。样地林木冠幅半径。单位：m。
#' @param Neighbor （必须）字符串。相邻木的选择方法，2种："number"（固定数量法，如k=4，6，8等）和"distance"（固定半径法）。
#' @param k （Neighbor="number"时必须）大于0小于样地林木数量的整数值。计算竞争指数的相邻木数量，默认为4。
#' @param Search_dist （Neighbor= "distance"时必须）数字向量。固定半径法的半径。单位：m。
#' @param Shape （必须）字符串。输入分析样地的形状，"rectangle"（矩形样地）或"circle"（圆形样地）。默认为"rectangle"。
#' @param Correct （必须）字符串。输入边缘矫正方法，"translation"、"buffer"或"single"。默认为"single"，逐株优化法，即以距对象木最近的第k株相邻木的距离进行判断；"translation"为平移法（8邻域缓冲区法）；"buffer"为固定距离的缓冲区法。
#' @param Origin （必须）数字向量，长度2。样地原点的坐标c(xmin, ymin)。默认为c(0, 0)。
#' @param Range_xy （Shape = "rectangle"时必填）数值向量，长度2，均大于0。矩形样地的横坐标X和纵坐标Y长度。单位：m。
#' @param Radius （Shape = "circle"时必须）数字，大于0。圆形样地的半径r。单位：m。
#' @param Buf_dist （correct = "buffer"和"translation"时必须）数字，大于0且小于短边的一半（Shape = "rectangle"时）或大于0且小于半径r（Shape = "circle"时）。缓冲区的宽度，默认为5m。单位：m。
#' @param Bind （必须）逻辑向量。是否将输出的空间结构指标与Data进行匹配合并。默认为FALSE，不进行合并。当Bind = TRUE时，将竞争指数按照样地号和林木标签与Data进行匹配后合并输出。

#' @return 1）正常输出结果为数据框。
#' @return 如果Bind=FLASE，结果数据框仅包含样地号、林木标签和竞争指数；若Bind=TRUE且Data不为空，结果数据框包含Data的所有林木信息和竞争指数。
#' @return 2）输出变量名：
#' @return Plot：样地号
#' @return Tag：林木标签
#' @return in_out：一个记录林木是否处于核心区的变量，包括"in"和"out"两类
#' @return CI15：与距离有关竞争指数CI15
#' @return CI15_intra：CI15的种内竞争
#' @return CI15_inter：CI15的种间竞争
#' @return CI17：与距离有关竞争指数CI17
#' @return CI17_intra：CI17的种内竞争
#' @return CI17_inter：CI17的种间竞争
#' @return CI18：与距离有关竞争指数CI18
#' @return CI18_intra：CI18的种内竞争
#' @return CI18_inter：CI18的种间竞争,
#' @return N_intra：相同树种的竞争木数量
#' @return N_inter：不同树种的竞争木数量

#' @name CI.dist.base

#' @references Alemdağ, I.S. (1978). Evaluation of some competition indexes for the prediction of diameter increment in planted white spruce. Information Report Forest Management Institute (Canada). FMR-X-108.
#' @references Hegyi, F. (1974). A simulation model for managing jack-pine stands. In International Union of Forestry Research Organizations, Proceedings of meeting in 1973, Growth Models for Tree and Stand Simulations. Edited by J. Fries. Royal College of Forestry, Stockholm, Stockholom, Sweden. Research Note No. 30, pp. 74-90.
#' @references Martin, G.L., Ek, A.R. (1984). A comparison of competition measures and growth models for predicting plantation red pine diameter and height growth. Forset Science, 30, 731-743.
#' @references Pommerening A, Stoyan D. Edge-correction needs in estimating indices of spatial forest structure[J]. Canadian Journal of Forest Research, 2006, 36(7): 1723-1739.

#' @export CI.dist.base

#' @import deldir
#' @import dplyr
#' @import plyr
#' @import stats

#' @examples ## 加载内置数据
#' @examples data(ForestStatTool)
#' @examples ##bigplot是一个100m×100m的矩形样地的林木信息，使用Data参数直接输入林木信息，输入的数据框中需要至少包括Plot，Tag，X，Y，D等变量名，SP可选。
#' @examples b1 <- CI.dist.base (Data = bigplot, Shape = "rectangle", Neighbor = "number", k = 4,
#' @examples                     Correct = "translation", Origin = c(0,0), Range_xy = c(100,100))
#' @examples head(b1)

#' @examples #cdf是一个半径13.82m的圆形样地的林木信息，使用Data参数直接输入林木信息，输入的数据框中需要至少包括Plot，Tag，X，Y，D等变量名，SP可选。
#' @examples b2 <- CI.dist.base (Data =cdf, Shape = "circle", Neighbor = "number", k = 6,
#' @examples                    Correct = "buffer", Origin = c(0,0), Radius = 13.6, Buf_dist = 2)
#' @examples head(b2)

CI.dist.base <- function(Data = NULL, Plot, Tag, X, Y, D, SP = NULL, H = NULL, CR = NULL,
                           Neighbor = "number", k = 4, Search_dist = NULL,
                           Shape = "rectangle", Correct = "single",
                           Origin = c(0,0), Range_xy = NULL, Radius = NULL, Buf_dist = 5, Bind = FALSE){
  if(!require("deldir")){
    install.packages("deldir")
    library(deldir)
  }
  if(!require("plyr")){
    install.packages("plyr")
    library(plyr)
  }
  if(!require("dplyr")){
    install.packages("dplyr")
    library(dplyr)
  }

# Shape--样地形状，2种：rectangle（方形）和circle（圆形）；
# Correct--边缘矫正方法，3种，single（逐株判断法）, buffer（缓冲区法） 和 translation（平移法），用“in”和“out”来区分
# x0, y0--对象木的坐标； Near--对象木的相邻木的数据框,并按照与对象木的距离升序排列，行数为Nmax
# Origin--样地原点坐标，长度为2的向量，如c(0,0)
# Range_xy-- 方形样地的x,y坐标边长，单位：米； Radius--圆形样地的半径，单位：米
# Buf_dist--缓冲区距离，单位：米
# Neighbor--竞争木选择方法，2种：number（固定数量法，如k=4，6，8等）和distance（固定半径法，由参数Search_dist指定）
# Search_dist--固定半径法选择竞争木时要指定的半径
Correction.f <- function(Shape=Shape, Correct=Correct, x0, y0, Near=NULL, Neighbor=Neighbor, Search_dist=Search_dist,
                       Origin=Origin,  RangeX=RangeX, RangeY=RangeY,
                       Buf_dist=Buf_dist, Radius=Radius) {
  # 如果是固定数量法，3种边缘矫正方式均可用
  if(Neighbor=="number"){
    #判断是否选择逐株边缘矫正方式(比较距第Nmax株相邻木的距离和距边界距离的大小)
    if(Shape == "rectangle"){
      distxy <- c(x0-Origin[1],y0-Origin[2],RangeX-(x0-Origin[1]),RangeY-(y0-Origin[2]))#对象点到四条边的最短距离
      if(Correct == "single"){
        Nmax=nrow(Near)
        if(Nmax == 0) {
          in_out = NA
        }else{
          x.max <- Near$X[Nmax] - x0
          y.max <- Near$Y[Nmax] - y0
          distmax <- sqrt(x.max^2 + y.max^2)#第Nmax株相邻木到对象点的距离
          in_out <-  ifelse(distmax <= min(distxy), "in", "out")
        }
      }
      if(Correct == "buffer"){
        in_out <-  ifelse(min(distxy) >= Buf_dist, "in", "out")
      }
      if(Correct == "translation"){
        in_out <- ifelse(distxy[1] >= 0 & distxy[2] >= 0 & distxy[3] > 0 & distxy[4] > 0, "in", "out")#点在左、下边保留，在右、上边剔除
      }
    }
    if(Shape == "circle"){
      xx <- Origin[1] - x0
      yy <- Origin[2] - y0
      distxy <- Radius-sqrt(xx^2 + yy^2)#对象点到圆边的最短距离
      if(Correct == "single"){
        Nmax=nrow(Near)
        if(Nmax == 0) {
          in_out = NA
        }else{
          x.max <- Near$X[Nmax] - x0
          y.max <- Near$Y[Nmax] - y0
          distmax <- sqrt(x.max^2 + y.max^2)#第Nmax株相邻木到对象点的距离
          in_out <-  ifelse(distmax <= distxy, "in", "out")
        }
      }
      if(Correct == "buffer"){
        in_out <-  ifelse(distxy >= Buf_dist, "in", "out")
      }
    }
  }
  # 如果是指定半径法，使用缓冲区法和平移法进行边缘矫正，其中 缓冲区法的距离 Buf_dist 直接继承参数Search_dist
  if(Neighbor=="distance"){
    Buf_dist=Search_dist
    #判断是否选择逐株边缘矫正方式(比较距第Nmax株相邻木的距离和距边界距离的大小)
    if(Shape == "rectangle"){
      distxy <- c(x0-Origin[1],y0-Origin[2],RangeX-(x0-Origin[1]),RangeY-(y0-Origin[2]))#对象点到四条边的最短距离
      if(Correct == "buffer"){
        in_out <-  ifelse(min(distxy) >= Buf_dist, "in", "out")
      }
      if(Correct == "translation"){
        in_out <- ifelse(distxy[1] >= 0 & distxy[2] >= 0 & distxy[3] > 0 & distxy[4] > 0, "in", "out")#点在左、下边保留，在右、上边剔除
      }
    }
    if(Shape == "circle"){
      xx <- Origin[1] - x0
      yy <- Origin[2] - y0
      distxy <- Radius-sqrt(xx^2 + yy^2)#对象点到圆边的最短距离
      if(Correct == "buffer"){
        in_out <-  ifelse(distxy >= Buf_dist, "in", "out")
      }
    }
  }
  in_out
}

###选取k株相邻木
# x0, y0--对象木的坐标
# datax--对象木所在的整个样地数据
# k--相邻木株数
NearN.f <- function(x0, y0, datax, k=k){
  x1 <- datax$X - x0
  y1 <- datax$Y - y0
  datax$dist <- sqrt(x1^2 + y1^2)#分别求算各个林木到对象木的距离
  # NearN <-datax[order(datax$dist),]
  NearN <- dplyr::arrange(datax, dist)#按照距离排序结果取自身和最近的k株相邻木的信息
  NearN <- NearN[1:(k+1),]
  return(NearN)
}
###选取距对象木距离小于半径Search_dist的相邻木
# x0, y0--对象木的坐标
# datax--对象木所在的整个样地数据
# Search_dist--固定半径法选择竞争木时要指定的半径
StepRad.f<-function(x0,y0,datax,Search_dist=Search_dist){
  x1<-datax$X-x0
  y1<-datax$Y-y0
  datax$dist<-sqrt(x1^2+y1^2)#分别求算各个林木到对象木的距离
  NearRad<-datax[datax$dist < Search_dist,]#取自身和半径内的所有相邻木的信息
  NearRad <- dplyr::arrange(NearRad,dist)#按照距离排序
  return(NearRad)
}

###平移法加入8邻域缓冲区
Translation.f <- function(plot0, RangeX=RangeX, RangeY=RangeY){
  ##中心为研究样地a1
  a1 <- plot0
  ##建立8邻域缓冲区
  for(m in 2:9){
    assign(paste0('a',m),plot0)
  }
  ##此处下标为X坐标列
  a2$X <- a2$X-RangeX
  a3$X <- a3$X-RangeX
  a4$X <- a4$X-RangeX
  a7$X <- a7$X+RangeX
  a8$X <- a8$X+RangeX
  a9$X <- a9$X+RangeX
  ##此处下标为Y坐标列
  a2$Y <- a2$Y-RangeY
  a5$Y <- a5$Y-RangeY
  a7$Y <- a7$Y-RangeY
  a4$Y <- a4$Y+RangeY
  a6$Y <- a6$Y+RangeY
  a9$Y <- a9$Y+RangeY
  ##组合为大样地
  bigplot <- rbind(a1,a2,a3,a4,a5,a6,a7,a8,a9)
  #判断是否有重复点
  dupxy <- duplicatedxy(bigplot$X,bigplot$Y)
  if(T %in% dupxy){
    bigplot$dupxy <- dupxy
    bigplot <- subset(bigplot,dupxy=="FALSE")
  }
  bigplot
}

#---------Hegyi指标计算
# x0, y0--对象木的坐标
# D0--对象木的胸径
# sp0--对象木的树种
# competition--对象木和竞争木组成的数据框
Hegyi.f<-function(x0,y0,D0,sp0,competition){
  N<-nrow(competition)
  if("SP" %in% colnames(competition)){
    if(N<=1){
      Hegyi<-rep(0,11)
    } else {
      # x0 <- competition$X[1]
      # y0 <- competition$Y[1]
      # D0<-competition$D[1]
      # sp0<-competition$SP[1]
      competition <- dplyr::mutate(competition,X0=x0,Y0=y0,
                            dist=sqrt((X-X0)^2+(Y-Y0)^2))
      Near<-competition[2:N,]
      Near_intra<-Near[Near$SP==sp0,]
      Near_inter<-Near[Near$SP!=sp0,]
      M1<-nrow(Near_intra)
      M2<-nrow(Near_inter)
      if(M1 == 0){
        CI15.1<-0
        CI17.1<-0
        CI18.1<-0
      } else {
        CI15.1<-sum(Near_intra$D/Near_intra$dist)/D0
        CI17.1<-pi*sum((Near_intra$dist*D0/(D0+Near_intra$D))^2*Near_intra$D/Near_intra$dist/sum(Near_intra$D/Near_intra$dist))
        CI18.1<-sum(Near_intra$D/D0*exp(Near_intra$dist*16/(D0+Near_intra$D)))
      }
      if(M2 == 0){
        CI15.2<-0
        CI17.2<-0
        CI18.2<-0
      } else {
        CI15.2<-sum(Near_inter$D/Near_inter$dist)/D0
        CI17.2<-pi*sum((Near_inter$dist*D0/(D0+Near_inter$D))^2*Near_inter$D/Near_inter$dist/sum(Near_inter$D/Near_inter$dist))
        CI18.2<-sum(Near_inter$D/D0*exp(Near_inter$dist*16/(D0+Near_inter$D)))
      }
      CI15<-sum(Near$D/Near$dist)/D0
      CI17<-pi*sum((Near$dist*D0/(D0+Near$D))^2*Near$D/Near$dist/sum(Near$D/Near$dist))
      CI18<-sum(Near$D/D0*exp(Near$dist*16/(D0+Near$D)))
      Hegyi<-c("CI15"=CI15,
               "CI15_intra"=CI15.1,
               "CI15_inter"=CI15.2,
               "CI17"=CI17,
               "CI17_intra"=CI17.1,
               "CI17_inter"=CI17.2,
               "CI18"=CI18,
               "CI18_intra"=CI18.1,
               "CI18_inter"=CI18.2,
               "N_intra"=M1,
               "N_inter"=M2)
    }
  } else {
    if(N<=1){
      Hegyi<-rep(0,3)
    } else {
      x0 <- competition$X[1]
      y0 <- competition$Y[1]
      D0<-competition$D[1]

      competition <- dplyr::mutate(competition,X0=x0,Y0=y0,
                            dist=sqrt((X-X0)^2+(Y-Y0)^2))
      Near<-competition[2:N,]
      if(nrow(Near) == 0){
        CI15<-0
        CI17<-0
        CI18<-0
      } else {
        CI15<-sum(Near$D/Near$dist)/D0
        CI17<-pi*sum((Near$dist*D0/(D0+Near$D))^2*Near$D/Near$dist/sum(Near$D/Near$dist))
        CI18<-sum(Near$D/D0*exp(Near$dist*16/(D0+Near$D)))
      }
      Hegyi<-c("CI15"=CI15,"CI17"=CI17,"CI18"=CI18)
    }
  }

  return(Hegyi)
}


  #参数合理性判断
  if(Shape %in% c("rectangle", "circle")){
    Shape=Shape
  }else{
    stop("Your choice for 'Shape' is not available, please choose 'rectangle' or 'circle'.")
  }
  if(Correct %in% c("single", "buffer", "translation")){
    Correct=Correct
  }else{
    stop("Your choice for 'Correct' is not available, please choose 'single', 'buffer' or 'translation'.")
  }
  if(Shape == "rectangle"){
    if(is.null(Range_xy)){
      stop("The rectangle plot needs x and y coordinate ranges (Range_xy).")
    }else{
      if(length(Range_xy)==2){
        if(!is.numeric(Range_xy)){
          stop("'Range_xy' must be numeric.")
        }
        RangeX <- Range_xy[1]
        RangeY <- Range_xy[2]
      }else{
        stop("The length of Range_xy should be 2.")
      }
      if(RangeX <= 0 || RangeY <= 0){
        stop("The coordinate ranges of the rectangle plot (Range_xy) must be greater than 0.")
      }else{
        if(Correct == "buffer"){
          if(Buf_dist <= 0 || Buf_dist >= min(RangeX, RangeY)/2){
            stop("Please set a reasonable Buf_dist which need greater than 0 and smaller than the minimum value between one-half of the coordinate ranges (Range_xy).")
          }
        }
        if(Correct == "translation"){
          if(Buf_dist <= 0 || Buf_dist > min(RangeX, RangeY)){
            stop("Please set a reasonable Buf_dist which need greater than 0 and smaller than the minimum value between the x and y coordinate ranges (Range_xy).")
          }
        }
      }
    }
  }
  if(Shape == "circle"){
    if(is.null(Radius)){
      stop("The circular plot needs a 'Radius'.")
    }else{
      if(!is.numeric(Radius)){
        stop("'Radius' must be numeric.")
      }
      if(Radius <= 0){
        stop("The 'Radius' of the circular plot must be greater than 0.")
      }else{
        if(Correct == "buffer"){
          if(Buf_dist <= 0 || Buf_dist >= Radius){
            stop("Please set a reasonable Buf_dist which need smaller than the 'Radius'.")
          }
        }
        if(Correct == "translation"){
          stop("The 'translation' method in parameter of 'Correct'is not suitable for cirlular plot, \n please try 'single' or 'buffer' for the parameter of 'Correct'.")
        }
      }
    }
  }
  if(Neighbor == "distance"){
    if(Correct == "single"){
      stop("The 'single' method in parameter of 'Correct'is not suitable for 'distance' method in parameter of 'Neighbor', \n please try 'buffer' or 'translation' for the parameter of 'Correct'.")
    }
  }


  #数据输入可选择数据框，也可以是向量
  #Plot, Tag, X, Y, D为必须输入指标
  if(is.null(Data)){
    data0 <- data.frame("Plot" = Plot, "Tag" = Tag, "X" = X, "Y" = Y, "D" = D)
    if(!is.null(H)){
      data0$H <- H
    }
    if(!is.null(CR)){
      data0$CR <- CR
    }
    if(!is.null(SP)){
      data0$SP <- SP
    }
    if(Bind==TRUE){
      stop("Missing 'Data', 'Bind = TRUE' is not applicable.")
    }
  } else { # 用Data输入时的检查
    if(!"Plot"%in%colnames(Data)){
      stop("Missing column name of 'Plot' in Data.")
    }
    if(!"Tag"%in%colnames(Data)){
      stop("Missing column name of 'Tag' in Data.")
    }
    if(!"X"%in%colnames(Data)){
      stop("Missing column name of 'X' in Data.")
    }
    if(!"Y"%in%colnames(Data)){
      stop("Missing column name of 'Y' in Data.")
    }
    if(!"D"%in%colnames(Data)){
      stop("Missing column name of 'D' in Data.")
    }
    if("SP"%in%colnames(Data)){
      data0 <- select(Data, Plot, Tag, X, Y, D, SP)
    }else{
      data0 <- select(Data, Plot, Tag, X, Y, D)
    }
  }

  # 统计NA的比例，何潇-2022-12-8
  Na.f=function(x){ round(sum(is.na(x))/length(x)*100, 2) }
  #主要变量合理性检测
  # 修改为变量有缺失值为打印，而不是停止，下面有删除NA的操作。何潇2022-12-8
  # 如果参数的类型不对，如D不是数值型，则停止计算。何潇2022-12-8
  if(sum(is.na(data0$Plot))!=0){
    print(paste0(Na.f(data0$Plot), "% Missing value (or NA) in 'Plot'."))
  }
  if(sum(is.na(data0$Tag))!=0){
    print(paste0(Na.f(data0$Tag), "% Missing value (or NA) in 'Tag'."))
  }
  if(sum(is.na(data0$SP))!=0){
    print(paste0(Na.f(data0$SP), "% Missing value (or NA) in 'SP'."))
  }
  ##
  if(!is.numeric(data0$X)){
    stop("'X' must be numeric.")
  }else if(sum(is.na(data0$X))!=0){
    print(paste0(Na.f(data0$X), "% Missing value (or NA) in 'X'."))
  }
  if(!is.numeric(data0$Y)){
    stop("'Y' must be numeric.")
  }else if(sum(is.na(data0$Y))!=0){
    print(paste0(Na.f(data0$Y), "% Missing value (or NA) in 'Y'."))
  }
  if(!is.numeric(data0$D)){
    stop("'D' must be numeric.")
  }else if(sum(is.na(data0$D))!=0){
    print(paste0(Na.f(data0$D), "% Missing value (or NA) in 'D'."))
  }
  # # H 和CR只有输入了才检查
  # if(!is.null(data0$H)){
  #   if(!is.numeric(data0$H)){
  #     stop("'H' must be numeric.")
  #   }else if(T %in% is.na(data0$H)){
  #     print(paste0(Na.f(data0$H), "% Missing value (or NA) in 'H'."))
  #   }
  # }
  # if(!is.null(data0$CR)){
  #   if(!is.numeric(data0$CR)){
  #     stop("'CR' must be numeric.")
  #   }else if(T %in% is.na(data0$CR)){
  #     print(paste0(Na.f(data0$CR), "% Missing value (or NA) in 'CR'."))
  #   }
  # }

  #删除NA值防止无法运算，何潇-2022-11-27
  data0 <- na.omit(data0)
  if(nrow(data0)==0){
    temp = apply(data0, 2, function(x) sum(!is.na(x))==0)
    tempName = names(which(temp == TRUE))
    if(length(tempName)==1){
      stop(paste0("The column of", tempName, "is NA."))
    }else{
      stop(paste0("The columns of", tempName, "are NA."))
    }
  }

  #统计样地个数
  N <- nlevels(factor(data0$Plot))
  # 增加进度条，何潇-2022-11-28
  print("Srart Calculateing")
  progress.bar <- plyr::create_progress_bar("text")  #plyr包中的create_progress_bar函数创建一个进度条
  progress.bar$init(N)   #设置任务数，几个样地

  data_all <- data.frame()
  for (i in 1:N) {
    datax <- dplyr::filter(data0,Plot==levels(factor(data0$Plot))[i])
    #重复标签检查
    if(T %in% duplicated(datax$Tag)){
      stop("There are tags duplicated in data, please use the function 'Tag_Divide' or 'Tag_Remove' to solve, and use the new column of 'newTag' to analyse.")
    }
    #重复坐标点检查
    if(T %in% duplicatedxy(datax$X,datax$Y)){
      stop("There are coordinates duplicated in data, please use the function 'Coord_Move' or 'Coord_Remove' to solve, and use the new columns of 'newX' and 'newY' to analyse.")
    }
    #检查参数k
    if(!k %in% c(1:(nrow(datax)-1))){
      stop("The parameter 'k' needs to be assigned an integer which greater than 0 and smaller than the number of trees in any plot.")
    }

    ##初始赋值
    Hegyi.data = data.frame()
    #平移法边缘矫正
    if(Correct == "translation"){
      datax <- Translation.f(datax, RangeX, RangeY)
      datax <- datax[datax$X > Origin[1] - Buf_dist & datax$X <= Origin[1] + RangeX + Buf_dist &
                       datax$Y > Origin[2] - Buf_dist & datax$Y <= Origin[2] + RangeY + Buf_dist,]
    }
    for (j in 1:nrow(datax)){
      x0 <- datax$X[j]
      y0 <- datax$Y[j]
      d0 <- datax$D[j]
      if(!is.null(data0$SP)){
        sp0 <- datax$SP[j]
        #提取相邻木信息
        if(Neighbor == "number"){
          Near0 <- NearN.f(x0, y0, datax, k)
        }else if(Neighbor == "distance"){
          Near0 <- StepRad.f(x0, y0, datax, Search_dist)
        }
        #调用hegyi.f计算
        Hegyi.data.j <- Hegyi.f(x0,y0,D0=d0,sp0=sp0,Near0)
        # 边缘矫正方法
        Hegyi.data.j <- as.data.frame(t(Hegyi.data.j))
        Hegyi.data.j$in_out <- Correction.f(Shape=Shape, Correct=Correct, x0=x0, y0=y0, Near=Near0, Neighbor=Neighbor, Search_dist=Search_dist,
                                          Origin=Origin,  RangeX=RangeX, RangeY=RangeY,
                                          Buf_dist=Buf_dist, Radius=Radius)
        Hegyi.data.j$Plot=levels(factor(data0$Plot))[i]
        Hegyi.data.j$Tag=datax$Tag[j]
        Hegyi.data <- rbind(Hegyi.data, Hegyi.data.j)
      }else {
        #提取相邻木信息
        if(Neighbor == "number"){
          Near0 <- NearN.f(x0, y0, datax, k)
        }else if(Neighbor == "distance"){
          Near0 <- StepRad.f(x0, y0, datax, Search_dist)
        }
        #调用hegyi.f计算
        Hegyi.data.j <- Hegyi.f(x0,y0,D0=d0,competition=Near0)
        # 边缘矫正方法
        Hegyi.data.j <- as.data.frame(t(Hegyi.data.j))
        Hegyi.data.j$in_out <- Correction.f(Shape=Shape, Correct=Correct, x0=x0, y0=y0, Near=Near0, Neighbor=Neighbor, Search_dist=Search_dist,
                                          Origin=Origin,  RangeX=RangeX, RangeY=RangeY,
                                          Buf_dist=Buf_dist, Radius=Radius)
        Hegyi.data.j$Plot=levels(factor(data0$Plot))[i]
        Hegyi.data.j$Tag=datax$Tag[j]
        Hegyi.data <- rbind(Hegyi.data, Hegyi.data.j)
      }
    }
    data_all <- rbind(data_all, Hegyi.data)
    # 每计算完一个样地打印一次结果
    print(paste("Calculated" ,i, ": Plot =", levels(factor(data0$Plot))[i]))
    progress.bar$step() #输出进度条
  }
  if(Correct == "translation"){
    data_all <- dplyr::filter(data_all, in_out == "in")
    }
  #判断是否与原数据框合并
  if(Bind){
    data_all <- dplyr::left_join(Data, data_all, by = c("Plot","Tag"))
  }
  return(data_all)
}

