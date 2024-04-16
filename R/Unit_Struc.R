#' @title 空间结构指数计算

#' @author Chaofan Zhou: cfzhou2021@163.com
#' @author Xiao He: hexiao@ifrit.ac.cn
#' @author Guangshuang Duan: oliverdgs@163.com
#' @author Jingning Shi: shijingning@hebau.edu.cn
#' @author Xiangdong Lei: xdlei@ifrit.ac.cn

#' @usage Unit_Struc(Data = NULL, Plot, Tag, X, Y, D, SP, H = NULL, CR = NULL, k = 4, Shape = "rectangle", Correct = "single", Origin = c(0,0), Range_xy = NULL, Radius = NULL, Buf_dist = 5, Bind = FALSE)

#' @description 此函数实现了根据结构单元（由一株中心木和其k株最近相邻木构成）计算林木和林分空间结构的功能。
#' @description 函数对矩形样地和圆形样地均适用（参数Shape），并内置3种边缘矫正方法（参数Correct），可同时分析多个样地的空间结构。

#' @details 1）输入数据可以选择两种方法：第一种是按照参数Plot, Tag, X, Y, D, SP, H, CR逐个输入其中H和CR可选择不输入；第二种是首先修改数据的列名，包括Plot, Tag, X, Y, D, SP, H, CR，其中H和CR可不包含在数据中，然后通过Data参数一次性输入。
#' @details 2）函数内置运算前检查功能，包括判断各个输入参数是否满足计算要求，是否存在重复坐标和标签。
#' @details 3）角尺度W_Upa仅在k大于等于3时计算。
#' @details 4）Buf_dist在correct = "translation"时的作用是确定原样地向外扩展的宽度，从而提高运算效率。

#' @param Data （可选）数据框。样地的所有林木属性信息，当Data输入时，Data中必须包含名称为Plot, Tag, X, Y, D, SP。在Bind = TRUE 时，Data必须输入。
#' @param Plot （必须）向量，建议为字符型。样地号。可以是数值型，也可以是字符型。
#' @param Tag （必须）向量。林木标签（编号）。可以是数值型，也可以是字符型。
#' @param X （必须）数值向量。样地林木横坐标。单位：m。
#' @param Y （必须）数值向量。样地林木纵坐标。单位：m。
#' @param D （必须）数值向量。样地林木胸径。单位：cm。
#' @param SP （必须）向量，可以是数值型，也可以是字符型。树种。
#' @param H （可选）数值向量。样地林木树高。单位：m。
#' @param CR （可选）数值向量。样地林木冠幅半径。单位：m。
#' @param k （必须）大于0小于样地林木数量的整数值。计算空间结构指标的相邻木数量，默认为4，即参照惠刚盈（2019）的相关研究结论，4株相邻木计算的空间结构值解释性更强。
#' @param Shape （必须）字符串。输入分析样地的形状，"rectangle"（矩形样地）或"circle"（圆形样地）。默认为"rectangle"。
#' @param Correct （必须）字符串。输入边缘矫正方法，"translation"、"buffer"或"single"。默认为"single"，即以距对象木最近的第k株相邻木的距离进行判断；"translation"为平移法（8邻域缓冲区法）；"buffer"为固定距离的缓冲区法。
#' @param Origin （必须）数值向量，长度2。样地原点的坐标c(xmin, ymin)。默认为c(0, 0)。
#' @param Range_xy （Shape = "rectangle"时必填）数值向量，长度2，均大于0。矩形样地的横坐标X和纵坐标Y长度。单位：m。
#' @param Radius （Shape = "circle"时必填）数值，大于0。圆形样地的半径r。单位：m。
#' @param Buf_dist （correct = "buffer"和"translation"时需要）数值，大于0且小于短边的一半（Shape = "rectangle"时）或大于0且小于半径r（Shape = "circle"时）。缓冲区的宽度，默认为5m。单位：m。
#' @param Bind （必须）逻辑向量。是否将输出的空间结构指标与Data进行匹配合并。默认为FALSE，不进行合并。当Bind = TRUE时，将空间结构指标按照样地号和林木标签与Data进行匹配后合并输出。


#' @return 1）正常输出结果为包含tree_value数据框和stand_value数据框的一个列表。
#' @return 如果Bind=FLASE，tree_value数据框仅包含样地号、林木标签和空间结构指标；若Bind=TRUE且Data不为空，tree_value数据框包含Data的所有林木信息和空间结构指标，若Data为空，则输出结果与Bind=FLASE时一致
#' @return stand_value数据框包含样地号和空间结构指标样地均值。
#' @return 2）输出变量名：
#' @return tree_value数据框中：
#' @return Plot：样地号
#' @return Tag：林木标签
#' @return in_out：一个记录林木是否处于核心区的变量，包括"in"和"out"两类
#' @return W_Upa：角尺度
#' @return U_Udm：直径大小比数
#' @return U_Uht：树高大小比数
#' @return U_Ucr：冠幅大小比数
#' @return T_Udm：直径大小分化度
#' @return T_Uht：树高大小分化度
#' @return T_Ucr：冠幅大小分化度
#' @return M_Usp：混交度
#' @return C_Ucr：密集度
#' @return … ：其他Data中的输入变量
#' @return stand_value数据框中：
#' @return Plot：样地号
#' @return Avg_W_Upa：角尺度的样地平均值
#' @return Avg_U_Udm：直径大小比数的样地平均值
#' @return Avg_U_Uht：树高大小比数的样地平均值
#' @return Avg_U_Ucr：冠幅大小比数的样地平均值
#' @return Avg_T_Udm：直径大小分化度的样地平均值
#' @return Avg_T_Uht：树高大小分化度的样地平均值
#' @return Avg_T_Ucr：冠幅大小分化度的样地平均值
#' @return Avg_M_Usp：混交度的样地平均值
#' @return Avg_C_Ucr：密集度的样地平均值

#' @references Pommerening A, Stoyan D. Edge-correction needs in estimating indices of spatial forest structure[J]. Canadian Journal of Forest Research, 2006, 36(7): 1723-1739.
#' @references Hui G, Zhang G, Zhao Z, et al. Methods of forest structure research: A review[J]. Current Forestry Reports, 2019, 5(3): 142-154.
#' @references Gadow K, González J, Zhang C, et al. Sustaining Forest Ecosystems[M].Springer.
#'
#' @export Unit_Struc
#' @name Unit_Struc

#' @import deldir
#' @import dplyr
#' @import plyr
#' @import stats

#' @examples ## 加载内置数据
#' @examples data(ForestStatTool)
#' @examples #使用参数逐个输入林木信息
#' @examples a1 <- Unit_Struc (Plot = bigplot$Plot, Tag = bigplot$Tag, X = bigplot$X, Y = bigplot$Y, D = bigplot$D,
#' @examples  SP = bigplot$SP, H = bigplot$H, Shape = "rectangle", Correct = "translation", Range_xy = c(100,100))
#' @examples head(a1[['tree_value']])
#' @examples head(a1[['stand_value']])

#' @examples #使用Data直接输入林木信息，Data中需要至少包括Plot，Tag，X，Y，D，SP等变量名，H和CR可选。
#' @examples #缓冲区法进行边缘矫正
#' @examples a2 <- Unit_Struc (Data = bigplot, Shape = "rectangle", Correct = "buffer",
#' @examples                   Range_xy = c(100,100), Buf_dist = 2, Bind = TRUE)
#' @examples head(a2[['tree_value']])
#' @examples head(a2[['stand_value']])

#' @examples #逐株优化法进行边缘矫正
#' @examples a3 <- Unit_Struc (Data = bigplot, k = 6, Shape = "rectangle", Correct = "single",
#' @examples                  Range_xy = c(100,100), Buf_dist = 5, Bind = TRUE)
#' @examples head(a3[['tree_value']])
#' @examples head(a3[['stand_value']])

#' @examples # cdf是一个圆形样地的林木信息
#' @examples a4 <- Unit_Struc (Plot = cdf$Plot, Tag = cdf$Tag, X = cdf$X, Y = cdf$Y, D = cdf$D, SP = cdf$SP,
#' @examples                  H = cdf$H, k = 3, Shape = "circle", Correct = "single", Radius = 13.82,
#' @examples                   Bind = FALSE)
#' @examples head(a4[['tree_value']])
#' @examples head(a4[['stand_value']])

Unit_Struc <- function(Data = NULL, Plot, Tag, X, Y, D, SP, H = NULL, CR = NULL,k = 4, Shape = "rectangle",Correct = "single",
                         Origin = c(0,0), Range_xy = NULL, Radius = NULL, Buf_dist = 5,  Bind = FALSE){
  # if(!require("deldir")){
    # install.packages("deldir")
    # library(deldir)
  # }
  # if(!require("plyr")){
    # install.packages("plyr")
    # library(plyr)
  # }
  # if(!require("dplyr")){
    # install.packages("dplyr")
    # library(dplyr)
  # }

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
          stop("The 'translation' method in parameter of 'Correct'is not suitable for cirlular plot, please try 'single' or 'buffer' for the parameter of 'Correct'.")
        }
      }
    }
  }

  # if(Shape == "polygon"){
  #   if(is.null(Clipxy)){
  #     stop("The polygon plot needs a 'Clipxy'.")
  #   }else{
  #     if(Radius <= 0){
  #       stop("The 'Radius' of the circular plot must be greater than 0.")
  #     }else{
  #       if(Correct == "buffer"){
  #         if(Buf_dist <= 0 || Buf_dist >= Radius){
  #           stop("Please set a reasonable Buf_dist which need smaller than the 'Radius'.")
  #         }
  #       }
  #       if(Correct == "translation"){
  #         stop("The 'translation' method in parameter of 'Correct'is not suitable for cirlular plot, please try 'single' or 'buffer' for the parameter of 'Correct'.")
  #       }
  #     }
  #   }
  # }


  #============Functions================
  #最近N株相邻木
  NearN.f <- function(x0, y0, datax, k){
    x1 <- datax$X - x0
    y1 <- datax$Y - y0
    datax$dist <- sqrt(x1*x1 + y1*y1)#分别求算各个林木到对象木的距离
    NearN <- datax[order(datax$dist),]#按照距离排序结果取自身和最近的k株相邻木的信息进行空间结构计算
    NearN <- NearN[1:(k+1),]
    return(NearN)
  }
  ##角尺度W
  Angle.f <- function(x0, y0, xi, yi){
    deltax <- xi-x0
    deltay <- yi-y0
    if (deltay > 0 & deltax >= 0) {Angle <- atan(deltax/deltay)} else
      if (deltay == 0 & deltax >= 0) {Angle <- pi*0.5} else
        if (deltay <= 0 & deltax >= 0) {Angle <- pi+atan(deltax/deltay)} else
          if (deltay < 0 & deltax < 0) {Angle <- pi+atan(deltax/deltay)} else
            if (deltay == 0 & deltax < 0) {Angle <- pi*1.5} else
            {Angle <- pi*2 + atan(deltax/deltay)}
    Angle
  }
  AngleDiff.f <- function(Angle1, Angle2, k){
    AngleDiff <- Angle1-Angle2
    if (AngleDiff > pi) {AngleDiff <- 2*pi - AngleDiff}
    if (AngleDiff >= 2*pi/(k+1)){AngleDiff <- 0} else
    {AngleDiff <- 1}
  }
  w.f <- function(x0, y0, x, y, k){
    Angle <- c()
    for (m in 1:k){
      xi <- x[m];yi <- y[m]
      Angle[m] <- Angle.f(x0, y0, xi, yi)
    }
    Order0 <- order(Angle)
    count <- 0
    for (m in 1:(k-1)) {
      count <- count + AngleDiff.f(Angle[Order0[m+1]], Angle[Order0[m]], k)
    }
    count <- count + AngleDiff.f(Angle[Order0[k]], Angle[Order0[1]], k)
    w <- count/k
  }
  ##大小比数U
  u.f <- function (Z0, Near, k, compare){
    if(compare =="D"){
      Z <- Near$D
    }else if(compare =="H"){
      Z <- Near$H
    }else if(compare =="CR"){
      Z <- Near$CR
    }
    count <- 0
    for (m in 1:k){
      if (Z[m] >= Z0){
        count <- count + 1
      }
    }
    u <- count/k
  }
  ##混交度M
  m.f <- function (sp0, Near, k){
    spi <- Near$SP
    count <- 0
    for (m in 1:k){
      if (spi[m] != sp0){
        count <- count + 1
      }
    }
    m <- count/k
  }
  ##密集度C
  c.f <- function (x0, y0, c0, Near, k){
    count <- 0
    for (m in 1:k){
      if (sqrt((Near$X[m]-x0)^2 + (Near$Y[m]-y0)^2) < (Near$CR[m]+c0)){
        count <- count + 1
      }
    }
    c <- count/k
  }
  ##大小分化度T
  t.f <- function (Z0, Near, k, compare){
    if(compare =="D"){
      Z <- Near$D
    }else if(compare =="H"){
      Z <- Near$H
    }else if(compare =="CR"){
      Z <- Near$CR
    }
    count <- 0
    for (m in 1:k){
      count <- count + min(Z[m],Z0)/max(Z[m],Z0)
    }
    t <- 1-count/k
  }
  ##平移法加入8邻域缓冲区
  Translation.f <- function(plot0, RangeX, RangeY){
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
    dupxy <- deldir::duplicatedxy(bigplot$X,bigplot$Y)
    if(T %in% dupxy){
      bigplot$dupxy <- dupxy
      bigplot <- subset(bigplot,dupxy=="FALSE")
    }
    bigplot
  }

  #数据输入可选择数据框，也可以是向量
  #Plot, Tag, X, Y, D, SP为必须输入指标
  if(is.null(Data)){
    data0 <- data.frame("Plot" = Plot, "Tag" = Tag, "X" = X, "Y" = Y, "D" = D, "SP" = SP)
    if(!is.null(H)){
      data0$H <- H
    }
    if(!is.null(CR)){
      data0$CR <- CR
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
    if(!"SP"%in%colnames(Data)){
      stop("Missing column name of 'SP' in Data.")
    }
    if("H"%in%colnames(Data)){
      if("CR"%in%colnames(Data)){
        data0 <- dplyr::select(Data, Plot, Tag, X, Y, D, SP, H, CR)
      }else{
        data0 <- dplyr::select(Data, Plot, Tag, X, Y, D, SP, H)
      }
    }else{
      if("CR"%in%colnames(Data)){
        data0 <- dplyr::select(Data, Plot, Tag, X, Y, D, SP, CR)
      }else{
        data0 <- dplyr::select(Data, Plot, Tag, X, Y, D, SP)
      }
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
  # H 和CR只有输入了才检查
  if(!is.null(data0$H)){
    if(!is.numeric(data0$H)){
      stop("'H' must be numeric.")
    }else if(T %in% is.na(data0$H)){
      print(paste0(Na.f(data0$H), "% Missing value (or NA) in 'H'."))
    }
  }
  if(!is.null(data0$CR)){
    if(!is.numeric(data0$CR)){
      stop("'CR' must be numeric.")
    }else if(T %in% is.na(data0$CR)){
      print(paste0(Na.f(data0$CR), "% Missing value (or NA) in 'CR'."))
    }
  }

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
    # if(TRUE %in% duplicated(datax$Tag)) {} # 何潇修改判断语句写法,2024-4-10
    if( any(duplicated(datax$Tag)) ){
      stop("There are tags duplicated in data, please use the function 'Tag_Divide' or 'Tag_Remove' to solve, and use the new column of 'newTag' to analyse.")
    }
    #重复坐标点检查
    # if(T %in% deldir::duplicatedxy(datax$X,datax$Y)){} # 何潇修改判断语句写法,2024-4-10
    if( any(deldir::duplicatedxy(datax$X, datax$Y)) ){
      stop("There are coordinates duplicated in data, please use the function 'Coord_Move' or 'Coord_Remove' to solve, and use the new columns of 'newX' and 'newY' to analyse.")
    }
    #检查参数k
    if(!k %in% c(1:(nrow(datax)-1))){
      stop("The parameter 'k' needs to be assigned an integer which greater than 0 and smaller than the number of trees in any plot.")
    }

    ##初始赋值
    W_Upa <- c()
    U_Udm <- c()
    U_Uht <- c()
    U_Ucr <- c()
    M_Usp <- c()
    C_Ucr <- c()
    T_Udm <- c()
    T_Uht <- c()
    T_Ucr <- c()
    in_out <- c()

    #平移法边缘矫正
    if(Correct == "translation"){
      datax <- Translation.f(datax, RangeX, RangeY)
      datax <- datax[datax$X > Origin[1] - Buf_dist & datax$X <= Origin[1] + RangeX + Buf_dist &
                     datax$Y > Origin[2] - Buf_dist & datax$Y <= Origin[2] + RangeY + Buf_dist,]
    }
    #循环计算空间结构指标
    for (j in 1:nrow(datax)){
      x0 <- datax$X[j]
      y0 <- datax$Y[j]
      d0 <- datax$D[j]
      sp0 <- datax$SP[j]
      #提取相邻木信息
      Near0 <- NearN.f(x0, y0, datax, k)
      Near <- Near0[2:(k+1),]
      #必算指标
      U_Udm[j] <- u.f(Z0 = d0, Near, k, compare = "D")
      T_Udm[j] <- t.f(Z0 = d0, Near, k, compare = "D")
      M_Usp[j] <- m.f(sp0, Near, k)

      #判断指标可否计算
      if(k>2){
        W_Upa[j] <- w.f(x0, y0, Near$X, Near$Y, k)
      }
      if(!is.null(data0$H)){
        h0 <- datax$H[j]
        U_Uht[j] <- u.f(Z0 = h0, Near, k, compare = "H")
        T_Uht[j] <- t.f(Z0 = h0, Near, k, compare = "H")
      }
      if(!is.null(data0$CR)){
        cr0 <- datax$CR[j]
        C_Ucr[j] <- c.f(x0, y0, cr0, Near, k)
        U_Ucr[j] <- u.f(Z0 = cr0, Near, k, compare = "CR")
        T_Ucr[j] <- t.f(Z0 = cr0, Near, k, compare = "CR")
      }
      #判断是否选择逐株边缘矫正方式(比较距第k株相邻木的距离和距边界距离的大小)
      if(Shape == "rectangle"){
        distxy <- c(x0-Origin[1],y0-Origin[2],RangeX-(x0-Origin[1]),RangeY-(y0-Origin[2]))#对象点到四条边的最短距离
        if(Correct == "single"){
          xk <- Near$X[k] - x0
          yk <- Near$Y[k] - y0
          distk <- sqrt(xk*xk + yk*yk)#第k株相邻木到对象点的距离
          in_out[j] <-  ifelse(distk <= min(distxy), "in", "out")
        }
        if(Correct == "buffer"){
          in_out[j] <-  ifelse(min(distxy) >= Buf_dist, "in", "out")
        }
        if(Correct == "translation"){
          in_out[j] <- ifelse(distxy[1] >= 0 & distxy[2] >= 0 & distxy[3] > 0 & distxy[4] > 0, "in", "out")#点在左、下边保留，在右、上边剔除
        }
      }
      if(Shape == "circle"){
        xx <- Origin[1] - x0
        yy <- Origin[2] - y0
        distxy <- Radius-sqrt(xx*xx + yy*yy)#对象点到圆边的最短距离
        if(Correct == "single"){
          xk <- Near$X[k] - x0
          yk <- Near$Y[k] - y0
          distk <- sqrt(xk*xk + yk*yk)
          in_out[j] <-  ifelse(distk <= distxy, "in", "out")
        }
        if(Correct == "buffer"){
          in_out[j] <-  ifelse(distxy >= Buf_dist, "in", "out")
        }
      }
    }
    result <- data.frame("Plot" = datax$Plot, "Tag" = datax$Tag,  "in_out" = in_out,
                         'U_Udm' = U_Udm, 'T_Udm' = T_Udm, 'M_Usp' = M_Usp)

    #判断加入其他指标
    if(k>2){
      result$W_Upa <- W_Upa
    }
    if(!is.null(data0$H)){
      result$U_Uht <- U_Uht
      result$T_Uht <- T_Uht
    }
    if(!is.null(data0$CR)){
      result$C_Ucr <- C_Ucr
      result$U_Ucr <- U_Ucr
      result$T_Ucr <- T_Ucr
    }
    data_all<-rbind(data_all,result)
    #数据分列整理后合并
    chr_data <- data_all[,colnames(data_all)%in%c("Plot","Tag","in_out")]
    int_data <- data_all[,!colnames(data_all)%in%c("Plot","Tag","in_out")]
    int_data <- round(int_data,2)
    data_all <- cbind(chr_data,int_data)
    # 每计算完一个样地打印一次结果
    print(paste("Calculated" ,i, ": Plot =", levels(factor(data0$Plot))[i]))
    progress.bar$step() #输出进度条
  }
  # 选择边缘矫正的方式
  if(Correct == "translation"){
    data_all <- dplyr::filter(data_all,in_out=="in")
  }
  if(Correct == "buffer"){
    print(paste0("The buffer distance is ",Buf_dist,"m."))
  }
  tree_in <- dplyr::filter(data_all,in_out=="in")

  #计算均值
  stand_value <- data.frame()
  for (i in 1:N) {
    plotx <- dplyr::filter(tree_in,Plot==levels(factor(data0$Plot))[i])
    plotx <- plotx[,!colnames(plotx)%in%c("Plot","Tag","in_out")]
    mean0 <- colMeans(plotx)
    mean0 <- data.frame(t(mean0))
    colnames(mean0) <- paste0("Avg_",colnames(mean0))
    mean0$Plot <- levels(factor(data0$Plot))[i]
    stand_value <- rbind(stand_value,mean0)
  }

  #判断是否与原数据框合并
  if(Bind){
    data_all <- dplyr::left_join(Data, data_all, by = c("Plot","Tag"))
  }

  output <- list("tree_value" = data_all, "stand_value" = stand_value)
  return(output)
}
