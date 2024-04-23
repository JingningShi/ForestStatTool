#' @title Stand factors. 林分因子

#' @author Xiao He: hexiaonuist@163.com
#' @author Chaofan Zhou: cfzhou2021@163.com
#' @author Guangshuang Duan: oliverdgs@163.com
#' @author Jingning Shi: shijingning@hebau.edu.cn
#' @author Xiangdong Lei: xdlei@ifrit.ac.cn

#' @usage StandFactor(D, Dmin = 5, Plot, S, TreeType = NULL, TreeClass = NULL, SP= NULL, H = NULL,
#' @usage D0 = 20, beta = 1.605, V = NULL, Biomass = NULL, Carbon = NULL)

#' @description 此函数实现了输入单木测树因子经统计输出常用林分（样地）因子.
#'
#' @details 1）只有当样地内每一株样木都有树高时才计算优势高;
#' @details 2）当Biomass或Carbon输入为一个向量时，输出结果变量为Bio或Car，如果输入N列数据框，则输出结果变量为Bio_1~ Bio_N或Car_1~ Car_N。

#' @param D （必须）数值向量。样地林木胸径。单位：cm。
#' @param Dmin 数值。起测胸径，默认为5cm。
#' @param Plot （必须）向量，可以是数值型，也可以是字符型。样地号。
#' @param S （必须）数值。样地面积。单位：ha。
#' @param SP （可选）向量，可以是数值型，也可以是字符型。树种，当SP为空或者输入的SP向量只有一个水平时，不会计算树种组成式和分树种的林分因子。
#' @param TreeType （可选）数值向量。检尺类型。12和14表示枯死木，13表示采伐木，其余值均认为是活立木，可输入全部检尺类型，程序将自动计算活立木（检尺类型不等于12、13、14和50）、采伐木（检尺类型为14或50）和枯死木（检尺类型13和15）的林分因子。默认为空，即输入的数据均为活立木。
#' @param TreeClass （可选）向量，可以是数值型，也可以是字符型。用于分组计算不同类型的林分因子，比如将样地内上层林木和下层林木的林分因子。当输入的TreeClass向量只有一个水平时，则不会计算。
#' @param H （可选）数值向量。单木树高。单位：m。用于计算林分平均高和优势高。林分平均高是样地内所有含树高值的单木树高平均值；林分优势高包括两种算法：1、样地内前n [n=round(S*100, 0)，S为样地面积，单位：ha] 株树高最大的优势木树高平均值；2、样地内前n [n=round(S*100, 0)，S为样地面积，单位：ha] 株胸径最大的优势木树高平均值。默认为空。
#' @param D0 数值。标准平均胸径。用于计算林分密度指数和可加林分密度指数，默认为20cm。
#' @param beta 数值。林分密度指数的斜率参数，默认为1.605。
#' @param V （可选）数值向量。单木材积。单位：m3，默认为空。
#' @param Biomass （可选）数值向量或数据框。单木生物量。单位：kg，默认为空。
#' @param Carbon （可选）数值向量或数据框。单木碳储量。单位：kg，默认为空。

#' @return 1）输出结果为一个列表，列表由数据框组成。如果输入了单木材积V，则列表中含有5个数据框：第1个数据框是常见林分因子；第2个数据框是断面积树种组成；第3个数据框是株数树种组成；第4个数据框是蓄积树种组成；第5个数据框是分树种的常见林分因子。如果单木材积V为空，则列表中含有4个数据框：第1个数据框是常见林分因子；第2个数据框是断面积树种组成；第3个数据框是株数树种组成；第4个数据框是分树种的常见林分因子。
#' @return 2）输出变量名：
#' @return Plot：样地号
#' @return Dg：活立木的林分断面积加权平均胸径
#' @return Ba：活立木的林分断面积
#' @return N；活立木的株数密度
#' @return SDI：活立木的林分密度指数
#' @return aSDI：活立木的可加林分密度指数
#' @return Dgd：优势胸径，样地内前n (n=S*100，S为样地面积，单位ha) 株胸径最大的优势木的林分断面积加权平均胸径
#' @return Hmean：林分平均高
#' @return Hmean_w_SP：按树种断面积加权的林分平均高
#' @return Ht：林分优势高1
#' @return Ht_Counts：林分优势高1的计算株数
#' @return Hd：林分优势高2
#' @return Hd_Counts：林分优势高2的计算株数
#' @return V：活立木的林分蓄积量
#' @return Bio 或Bio_X	：活立木的林分生物量
#' @return Car或Car_X：活立木的林分碳储量
#' @return Cut_Dg：采伐木的林分断面积加权平均胸径
#' @return Cut_Ba：采伐木的林分断面积
#' @return Cut_N：采伐木的株数密度
#' @return Cut_SDI：采伐木的林分密度指数
#' @return Cut_aSDI：采伐木的可加林分密度指数
#' @return Dead_Dg：枯死木的林分断面积加权平均胸径
#' @return Dead_Ba：枯死木的林分断面积
#' @return Dead_N：枯死木的株数密度
#' @return Dead_SDI：枯死木的林分密度指数
#' @return Dead_aSDI：枯死木的可加林分密度指数
#' @return Ba.p_XXX：XXX树种的断面积占比，XXX根据输入的单木树种SP决定
#' @return N.p_XXX：XXX树种的株数占比，XXX根据输入的单木树种SP决定
#' @return V.p_XXX：XXX树种的蓄积占比，XXX根据输入的单木树种SP决定
#' @return Composition：树种组成式，括号中的数字为占比

#' @export StandFactor
#' @name StandFactor

#' @importFrom dplyr group_by slice_max left_join mutate filter select arrange
#' @importFrom tidyr pivot_longer
#' @import stats

#' @examples ## 加载内置数据
#' @examples data(ForestStatTool)
#' @examples #bigplot是一个100m×100m的矩形样地的林木信息，首先利用Plot_Divide函数将bigplot划分为后的25个20m×20m的小样地的合并数据框。
#' @examples #除原始变量外，还包含的新列有：分割后的新样地名subplot，小样地的新X坐标subX和小样地的新Y坐标subY
#' @examples subplot1 <- Plot_Divide(Data = bigplot, Plot = bigplot$Plot, X = bigplot$X, Y = bigplot$Y,
#' @examples                        Num_xy = c(5,5), Range_xy = c(100,100))
#' @examples #计算不同小样地的林分因子
#' @examples #输出了4个数据框，分别是样地林分因子、断面积占比与组成、株数占比与组成
#' @examples a1 = StandFactor(Plot = subplot1$subplot, SP=subplot1$SP, S=0.04, D=subplot1$D)
#' @examples head(a1[[1]])
#' @examples head(a1[[2]])
#' @examples head(a1[[3]])
#' @examples head(a1[[4]])

#' @examples #cdf是一个半径13.82m的圆形样地的林木信息，计算圆形样地的林分因子
#' @examples #输出了5个数据框，分别是样地林分因子、断面积占比与组成、株数占比与组成、蓄积占比与组成
#' @examples a2 = StandFactor(Plot = cdf$Plot, SP=cdf$SP, S=0.06,
#' @examples                  D=cdf$D, H=cdf$H, V=cdf$V,
#' @examples                  Biomass=data.frame(cdf$BIO_A,cdf$BIO_ROOT),
#' @examples                  Carbon=data.frame(cdf$CAR_A, cdf$CAR_ROOT))

#' @examples head(a2[[1]])
#' @examples head(a2[[2]])
#' @examples head(a2[[3]])
#' @examples head(a2[[4]])
#' @examples head(a2[[5]])



StandFactor <- function(D, Dmin = 5, Plot, S, TreeType = NULL, TreeClass = NULL, SP = NULL, H = NULL,
                        D0 = 20, beta = 1.605, V = NULL, Biomass = NULL, Carbon = NULL) {
  # 检测dplyr包是否存在，只有不存在时才会安装
  if (!require("dplyr")) {
    install.packages("dplyr")
    library(dplyr)
  }
  # 检测tidyr包是否存在，只有不存在时才会安装
  if (!require("tidyr")) {
    install.packages("tidyr")
    library(tidyr)
  }

  # 林分因子基本计算函数
  Dg <- function(x) {
    Dg <- round(sqrt(sum(x^2) / length(x)), 1)
  }
  Ba <- function(x) {
    Ba <- round(sum((pi * x^2 / 40000)) / S, 2)
  }
  N <- function(x) {
    N <- round(length(x) / S, 0)
  }
  SDI <- function(x) {
    D <- Dg(x)
    N <- N(x)
    SDI <- round(N * (D / D0)^beta, 0)
    return(SDI)
  }
  aSDI <- function(x) {
    aSDI <- round(sum((x / D0)^beta) / S, 0)
  }
  tempFUN <- function(x) {
    Dg <- Dg(x)
    Ba <- Ba(x)
    N <- N(x)
    SDI <- SDI(x)
    aSDI <- aSDI(x)
    result <- data.frame("Dg" = Dg, "Ba" = Ba, "N" = N, "SDI" = SDI, "aSDI" = aSDI)
  }
  # 树种组成计算函数
  CompositionStat <- function(zhanbi, start = 6) {
    zhanbi_long <- pivot_longer(zhanbi,
      cols = colnames(zhanbi[, -1]),
      names_to = "key",
      values_to = "ratio"
    )
    zhanbi_long <- subset(zhanbi_long, zhanbi_long$ratio != 0)
    zhanbi_long <- zhanbi_long %>% arrange(Plot, desc(ratio))
    zhanbi_long$SP <- substr(zhanbi_long$key, start, nchar(zhanbi_long$key))
    comp <- c()
    comp.Plot <- c()
    for (i in 1:length(unique(zhanbi_long$Plot))) {
      temp <- subset(zhanbi_long, zhanbi_long$Plot == unique(zhanbi_long$Plot)[i])
      temp.comp <- ""
      temp.Plot <- ""
      if (nrow(temp) == 1) {
        comp[i] <- paste0(temp$ratio, "(", temp$SP, ")")
        comp.Plot[i] <- temp$Plot
      } else {
        for (j in 1:nrow(temp)) {
          if (temp$ratio[j] < 0.2) {
            temp.comp <- paste0(temp.comp, "-", "(", temp$SP[j], ")")
          } else if (temp$ratio[j] >= 0.2 & temp$ratio[j] < 0.5) {
            temp.comp <- paste0(temp.comp, "+", "(", temp$SP[j], ")")
          } else {
            temp.comp <- paste0(temp.comp, temp$ratio[j], "(", temp$SP[j], ")")
          }
        }
        comp[i] <- temp.comp
        comp.Plot[i] <- as.character(unique(temp$Plot))
      }
    }
    Composition <- data.frame("Plot" = comp.Plot, "Composition" = comp)
    return(Composition)
  }
  # 统计NA的比例，何潇-2022-12-8
  Na.f=function(x){ round(sum(is.na(x))/length(x)*100, 2) }
  Na.muti.f=function(data){
    Na.temp=apply(data, 2, Na.f)
    Na.temp
  }
  # 输入检查
  if(is.null(Plot)){
    stop("Missing column name of 'Plot'.")
  }else { Plot <- as.character(Plot) }
  if(!is.numeric(D)){
    stop("'D' must be numeric.")
  }else if(sum(is.na(D))!=0){
    cat(paste0(Na.f(D), "% Missing value (or NA) in 'D'.\n"))
  }
  ba <- pi * D^2 / 40000
  if (!is.null(TreeType)) {
    TT <- as.character(TreeType)
    if(sum(is.na(TT))!=0){
      cat(paste0(Na.f(TT), "% Missing value (or NA) in 'TreeType'.\n"))
    }
  } else {
    TT <- rep(11, length(D))
  }
  if (!is.null(SP)) {
    SP <- as.character(SP)
    if(sum(is.na(SP))!=0){
      cat(paste0(Na.f(SP), "% Missing value (or NA) in 'SP'.\n"))
    }
    if (length(unique(SP)) == 1) {
      cat("The 'SP' was only one level. So, the stand factors of the tree species will not be calculated!\n")
    }
  } else {
    SP <- rep("Pure", length(D))
  }

  # 检查H
  if(!is.null(H)){
    if(!is.numeric(H)){
      stop("'H' must be numeric.")
    }else if(sum(is.na(H))!=0){
      cat(paste0(Na.f(H), "% Missing value (or NA) in 'H'.\n"))
    }
  }
  # 检查V
  if(!is.null(V)){
    if(!is.numeric(V)){
      stop("'V' must be numeric.")
    }else if(sum(is.na(V))!=0){
      cat(paste0(Na.f(V), "% Missing value (or NA) in 'V'.\n"))
    }
  }

  # 检查Biomass和Carbon
  is.numeric.muti.f=function(data){
    is.numeric.temp=apply(data, 2, is.numeric)
    is.numeric.temp
  }
  # 检查Biomass
  if(!is.null(Biomass)){
    if (ncol(as.data.frame(Biomass)) != 1){
      if(sum(is.numeric.muti.f(Biomass))==0){
        stop("Each 'Biomass' must be numeric.")
      } else if(sum(Na.muti.f(Biomass))!=0){
        for (i.col in 1:ncol(as.data.frame(Biomass))) {
          cat(paste0(Na.muti.f(Biomass)[i.col], "% Missing value (or NA) in 'Biomass' col", i.col, ".\n"))
        }
      }
    } else{
      if(!is.numeric(Biomass)){
        stop("'Biomass' must be numeric.")
      }else if(sum(is.na(Biomass))!=0){
        cat(paste0(Na.f(Biomass), "% Missing value (or NA) in 'Biomass'.\n"))
      }
    }
  }
  # 检查Carbon
  if(!is.null(Carbon)){
    if (ncol(as.data.frame(Carbon)) != 1){
      if(sum(is.numeric.muti.f(Carbon))==0){
        stop("Each 'Carbon' must be numeric.")
      } else if(sum(Na.muti.f(Carbon))!=0){
        for (i.col in 1:ncol(as.data.frame(Carbon))) {
          cat(paste0(Na.muti.f(Carbon)[i.col], "% Missing value (or NA) in 'Carbon' col", i.col, ".\n"))
        }
      }
    } else{
      if(!is.numeric(Carbon)){
        stop("'Carbon' must be numeric.")
      }else if(sum(is.na(Carbon))!=0){
        cat(paste0(Na.f(Carbon), "% Missing value (or NA) in 'Carbon'.\n"))
      }
    }
  }

  # 数据组织
  data <- data.frame("Plot" = Plot, "SP" = SP, "TreeType" = TT, "d" = D, "ba" = ba)

  #  活立木大条件：!检尺类型%in%c(13,14,15,50)
  data.Alive <- subset(data, (data$d >= Dmin) & (!data$TreeType %in% c(13, 14, 15, 50)))

  if(nrow(data.Alive)==0){
    temp = apply(data.Alive, 2, function(x) sum(!is.na(x))==0)
    tempName = names(which(temp == TRUE))
    if(length(tempName)==1){
      stop(paste0("The column of", tempName, "is NA."))
    }else{
      stop(paste0("The columns of", tempName, "are NA."))
    }
  }

  Alive.StandFactor <- tapply(data.Alive$d, INDEX = data.Alive$Plot, FUN = tempFUN)
  Alive.StandFactor2 <- data.frame(matrix(unlist(Alive.StandFactor, use.names = T), ncol = 5, byrow = T), stringsAsFactors = F)
  colnames(Alive.StandFactor2) <- c("Dg", "Ba", "N", "SDI", "aSDI")

  # 活立木大的基本林分因子
  Alive.StandFactor3 <- data.frame("Plot" = rownames(Alive.StandFactor), Alive.StandFactor2)

  # 优势胸径，2022-1-12
  data.dominantDg <- data.Alive %>%
    group_by(Plot) %>%
    slice_max(n = round(S * 100, 0), order_by = d, with_ties = F)
  dominantDg <- tapply(data.dominantDg$d, INDEX = data.dominantDg$Plot, FUN = Dg)
  dominantDg <- data.frame("Plot" = rownames(dominantDg), "Dgd" = dominantDg)
  Alive.StandFactor3 <- left_join(Alive.StandFactor3, dominantDg, by = "Plot")

  # 树种组成，只有输入的 SP 不止一个树种时才计算，2021-10-16
  if (length(unique(data.Alive$SP)) != 1) {
    # 断面积比例----单独作为一个文件比较好，方便数据合并
    ba.sum <- tapply(data.Alive$ba, INDEX = list(data.Alive$Plot, data.Alive$SP), FUN = sum)
    ba.sum[is.na(ba.sum)] <- 0
    ba.zhanbi <- round(ba.sum / rowSums(ba.sum) * 100, 0)
    ba.zhanbi[is.na(ba.zhanbi)] <- 0
    ba.zhanbi[sapply(ba.zhanbi, is.infinite)] <- 0
    # 归一化
    for (i in 1:nrow(ba.zhanbi)) {
      zuchenghe <- sum(ba.zhanbi[i, ])
      if (zuchenghe != 100) {
        if (zuchenghe < 100) {
          ba.zhanbi[i, which.max(ba.zhanbi[i, ])] <- ba.zhanbi[i, which.max(ba.zhanbi[i, ])] + (100 - zuchenghe)
        } else {
          ba.zhanbi[i, which.max(ba.zhanbi[i, ])] <- ba.zhanbi[i, which.max(ba.zhanbi[i, ])] - (zuchenghe - 100)
        }
      }
    }
    # 重命名
    ba.names <- colnames(ba.zhanbi)
    ba.names <- paste0("Ba.p_", ba.names)
    colnames(ba.zhanbi) <- ba.names
    ba.zhanbi <- data.frame("Plot" = rownames(as.data.frame(ba.sum)), ba.zhanbi)
    # 断面积树种组成式
    ba.Composition <- CompositionStat(ba.zhanbi)
    ba.zhanbi <- left_join(ba.zhanbi, ba.Composition, by = "Plot")

    # 株数比例----单独作为一个文件比较好，方便数据合并
    counts.sp <- tapply(data.Alive$d, INDEX = list(data.Alive$Plot, data.Alive$SP), FUN = length)
    counts.sp[is.na(counts.sp)] <- 0
    n.zhanbi <- round(counts.sp / rowSums(counts.sp) * 100, 0)
    n.zhanbi[is.na(n.zhanbi)] <- 0
    n.zhanbi[sapply(n.zhanbi, is.infinite)] <- 0
    # 归一化
    for (i in 1:nrow(n.zhanbi)) {
      zuchenghe <- sum(n.zhanbi[i, ])
      if (zuchenghe != 100) {
        if (zuchenghe < 100) {
          n.zhanbi[i, which.max(n.zhanbi[i, ])] <- n.zhanbi[i, which.max(n.zhanbi[i, ])] + (100 - zuchenghe)
        } else {
          n.zhanbi[i, which.max(n.zhanbi[i, ])] <- n.zhanbi[i, which.max(n.zhanbi[i, ])] - (zuchenghe - 100)
        }
      }
    }
    # 重命名
    n.names <- colnames(n.zhanbi)
    n.names <- paste0("N.p_", n.names)
    colnames(n.zhanbi) <- n.names
    n.zhanbi <- data.frame("Plot" = rownames(as.data.frame(counts.sp)), n.zhanbi)
    # 株数树种组成式
    n.Composition <- CompositionStat(n.zhanbi, start = 5)
    n.zhanbi <- left_join(n.zhanbi, n.Composition, by = "Plot")
  }

  # 分树种计算林分因子
  if (length(unique(data.Alive$SP)) != 1) {
    # Dg
    Alive.SP.Dg <- tapply(data.Alive$d, INDEX = list(data.Alive$Plot, data.Alive$SP), FUN = Dg)
    Alive.SP.Dg[is.na(Alive.SP.Dg)] <- 0
    Alive.SP.Dg.names <- colnames(Alive.SP.Dg)
    Alive.SP.Dg.names <- paste0("Dg_", Alive.SP.Dg.names)
    colnames(Alive.SP.Dg) <- Alive.SP.Dg.names
    # Ba
    Alive.SP.Ba <- tapply(data.Alive$d, INDEX = list(data.Alive$Plot, data.Alive$SP), FUN = Ba)
    Alive.SP.Ba[is.na(Alive.SP.Ba)] <- 0
    Alive.SP.Ba.names <- colnames(Alive.SP.Ba)
    Alive.SP.Ba.names <- paste0("Ba_", Alive.SP.Ba.names)
    colnames(Alive.SP.Ba) <- Alive.SP.Ba.names
    # N
    Alive.SP.N <- tapply(data.Alive$d, INDEX = list(data.Alive$Plot, data.Alive$SP), FUN = N)
    Alive.SP.N[is.na(Alive.SP.N)] <- 0
    Alive.SP.N.names <- colnames(Alive.SP.N)
    Alive.SP.N.names <- paste0("N_", Alive.SP.N.names)
    colnames(Alive.SP.N) <- Alive.SP.N.names
    # SDI
    Alive.SP.SDI <- tapply(data.Alive$d, INDEX = list(data.Alive$Plot, data.Alive$SP), FUN = SDI)
    Alive.SP.SDI[is.na(Alive.SP.SDI)] <- 0
    Alive.SP.SDI.names <- colnames(Alive.SP.SDI)
    Alive.SP.SDI.names <- paste0("SDI_", Alive.SP.SDI.names)
    colnames(Alive.SP.SDI) <- Alive.SP.SDI.names
    # aSDI
    Alive.SP.aSDI <- tapply(data.Alive$d, INDEX = list(data.Alive$Plot, data.Alive$SP), FUN = aSDI)
    Alive.SP.aSDI[is.na(Alive.SP.aSDI)] <- 0
    Alive.SP.aSDI.names <- colnames(Alive.SP.aSDI)
    Alive.SP.aSDI.names <- paste0("aSDI_", Alive.SP.aSDI.names)
    colnames(Alive.SP.aSDI) <- Alive.SP.aSDI.names

    # 分树种的基本林分因子
    Alive.SP <- data.frame("Plot" = rownames(Alive.SP.Dg), Alive.SP.Dg, Alive.SP.Ba, Alive.SP.N, Alive.SP.SDI, Alive.SP.aSDI)

    # 分树种的蓄积 2021-10-16
    if (!is.null(V)) {
      V <- as.numeric(V)
      data <- data.frame("Plot" = Plot, "SP" = SP, "TreeType" = TT, "d" = D, "V" = V)
      Alive.SP.V <- subset(data, (data$d >= Dmin) & (!data$TreeType %in% c(13, 14, 15, 50)))
      Alive.SP.V <- tapply(Alive.SP.V$V / S, INDEX = list(data.Alive$Plot, data.Alive$SP), FUN = sum)
      Alive.SP.V[is.na(Alive.SP.V)] <- 0
      Alive.SP.V.names <- colnames(Alive.SP.V)
      Alive.SP.V.names <- paste0("V_", Alive.SP.V.names)
      colnames(Alive.SP.V) <- Alive.SP.V.names
      Alive.SP.V <- data.frame("Plot" = rownames(Alive.SP.V), Alive.SP.V)
      # 拼接
      Alive.SP <- left_join(Alive.SP, Alive.SP.V, by = "Plot")
    }
    # 分树种的生物量 2021-10-17
    if (!is.null(Biomass)) {
      Biomass <- as.data.frame(apply(Biomass, MARGIN = 2, FUN = as.numeric))
      data <- data.frame("Plot" = Plot, "SP" = SP, "TreeType" = TT, "d" = D, Biomass)
      data.SP.Bio <- subset(data, (data$d >= Dmin) & (!data$TreeType %in% c(13, 14, 15, 50)))
      if (ncol(as.data.frame(Biomass)) != 1) {
        SP.StandBio_1 <- tapply(data.SP.Bio[, 4 + 1] / 1000 / S, INDEX = list(data.SP.Bio$Plot, data.SP.Bio$SP), FUN = sum)
        SP.StandBio_1[is.na(SP.StandBio_1)] <- 0
        SP.StandBio_1.names <- colnames(SP.StandBio_1)
        SP.StandBio_1.names <- paste0("Bio_1_", SP.StandBio_1.names)
        colnames(SP.StandBio_1) <- SP.StandBio_1.names
        SP.StandBio_1 <- data.frame("Plot" = rownames(SP.StandBio_1), SP.StandBio_1)
        for (i in 2:ncol(Biomass)) {
          SP.StandBio_i <- tapply(data.SP.Bio[, 4 + i] / 1000 / S, INDEX = list(data.SP.Bio$Plot, data.SP.Bio$SP), FUN = sum)
          SP.StandBio_i[is.na(SP.StandBio_i)] <- 0
          SP.StandBio_i.names <- colnames(SP.StandBio_i)
          colnames(SP.StandBio_i) <- paste0("Bio_", i, "_", SP.StandBio_i.names)
          SP.StandBio_1 <- cbind(SP.StandBio_1, SP.StandBio_i)
        }
        Alive.SP <- left_join(Alive.SP, SP.StandBio_1, by = "Plot")
      } else {
        SP.StandBio_1 <- tapply(data.SP.Bio[, 4 + 1] / 1000 / S, INDEX = list(data.SP.Bio$Plot, data.SP.Bio$SP), FUN = sum)
        SP.StandBio_1[is.na(SP.StandBio_1)] <- 0
        SP.StandBio_1.names <- colnames(SP.StandBio_1)
        SP.StandBio_1.names <- paste0("Bio_", SP.StandBio_1.names)
        colnames(SP.StandBio_1) <- SP.StandBio_1.names
        SP.StandBio_1 <- data.frame("Plot" = rownames(SP.StandBio_1), SP.StandBio_1)
        Alive.SP <- left_join(Alive.SP, SP.StandBio_1, by = "Plot")
      }
    }
    # 分树种的碳储量 2021-10-17
    if (!is.null(Carbon)) {
      Carbon <- as.data.frame(apply(Carbon, MARGIN = 2, FUN = as.numeric))
      data <- data.frame("Plot" = Plot, "SP" = SP, "TreeType" = TT, "d" = D, Carbon)
      data.SP.Car <- subset(data, (data$d >= Dmin) & (!data$TreeType %in% c(13, 14, 15, 50)))
      if (ncol(as.data.frame(Carbon)) != 1) {
        SP.StandCar_1 <- tapply(data.SP.Car[, 4 + 1] / 1000 / S, INDEX = list(data.SP.Car$Plot, data.SP.Car$SP), FUN = sum)
        SP.StandCar_1[is.na(SP.StandCar_1)] <- 0
        SP.StandCar_1.names <- colnames(SP.StandCar_1)
        SP.StandCar_1.names <- paste0("Car_1_", SP.StandCar_1.names)
        colnames(SP.StandCar_1) <- SP.StandCar_1.names
        SP.StandCar_1 <- data.frame("Plot" = rownames(SP.StandCar_1), SP.StandCar_1)
        for (i in 2:ncol(Carbon)) {
          SP.StandCar_i <- tapply(data.SP.Car[, 4 + i] / 1000 / S, INDEX = list(data.SP.Car$Plot, data.SP.Car$SP), FUN = sum)
          SP.StandCar_i[is.na(SP.StandCar_i)] <- 0
          SP.StandCar_i.names <- colnames(SP.StandCar_i)
          colnames(SP.StandCar_i) <- paste0("Car_", i, "_", SP.StandCar_i.names)
          SP.StandCar_1 <- cbind(SP.StandCar_1, SP.StandCar_i)
        }
        Alive.SP <- left_join(Alive.SP, SP.StandCar_1, by = "Plot")
      } else {
        SP.StandCar_1 <- tapply(data.SP.Car[, 4 + 1] / 1000 / S, INDEX = list(data.SP.Car$Plot, data.SP.Car$SP), FUN = sum)
        SP.StandCar_1[is.na(SP.StandCar_1)] <- 0
        SP.StandCar_1.names <- colnames(SP.StandCar_1)
        SP.StandCar_1.names <- paste0("Car_", SP.StandCar_1.names)
        colnames(SP.StandCar_1) <- SP.StandCar_1.names
        SP.StandCar_1 <- data.frame("Plot" = rownames(SP.StandCar_1), SP.StandCar_1)
        Alive.SP <- left_join(Alive.SP, SP.StandCar_1, by = "Plot")
      }
    }
  }

  ## 林分高计算
  ## 数据组织
  if (!is.null(H)) {
    data.h <- data.frame("Plot" = Plot, "TreeType" = TT, "d" = D, "h" = H)
    data.Alive.h <- subset(data.h, data.h$d >= Dmin & data.h$h > 1.3 & !data.h$TreeType %in% c(13, 14, 15, 50))
    # 样地算术平均高
    AverageH <- tapply(data.Alive.h$h, INDEX = data.Alive.h$Plot, FUN = mean)
    AverageH[is.na(AverageH)] <- 0
    AverageH <- round(AverageH, 1)
    AverageH <- data.frame("Plot" = rownames(AverageH), "Hmean" = AverageH)
    Alive.StandFactor3 <- left_join(Alive.StandFactor3, AverageH, by = "Plot")

    # 样地按径阶加权平均高
    # data.Alive.h$dClass = cut(data.Alive.h$d, breaks = c(2*(4:27)-3,300), labels = c(2*(4:26)-2, '>51'), right = F)
    # x = tapply(data.Alive.h$ba, INDEX=list(data.Alive.h$Plot, data.Alive.h$dClass), FUN=sum)
    # x[is.na(x)] = 0
    # x = x/rowSums(x)
    # y = tapply(data.Alive.h$h, INDEX=list(data.Alive.h$Plot, data.Alive.h$dClass), FUN=mean)
    # y[is.na(y)] = 0
    # Hmean_w_dClass = rowSums(x*y)
    # Hmean_w_dClass = data.frame('Plot'=names(Hmean_w_dClass), 'Hmean_w_dClass'=Hmean_w_dClass)
    # Alive.StandFactor3 = left_join(Alive.StandFactor3, Hmean_w_dClass, by='Plot')
    # 样地按树种断面积加权平均高
    if (length(unique(data.Alive$SP)) != 1) {
      data.h1 <- data.frame("Plot" = Plot, "TreeType" = TT, "d" = D, "ba" = ba, "h" = H, "SP" = SP)
      data.Alive.h1 <- subset(data.h1, data.h$d >= Dmin & data.h$h > 1.3 & !data.h$TreeType %in% c(13, 14, 15, 50))
      x <- tapply(data.Alive.h1$ba, INDEX = list(data.Alive.h1$Plot, data.Alive.h1$SP), FUN = sum)
      x[is.na(x)] <- 0
      x <- x / rowSums(x)
      y <- tapply(data.Alive.h1$h, INDEX = list(data.Alive.h1$Plot, data.Alive.h1$SP), FUN = mean)
      y[is.na(y)] <- 0
      Hmean_w_SP <- round(rowSums(x * y), 1)
      Hmean_w_SP <- data.frame("Plot" = names(Hmean_w_SP), "Hmean_w_SP" = Hmean_w_SP)
      Alive.StandFactor3 <- left_join(Alive.StandFactor3, Hmean_w_SP, by = "Plot")
    }

    # 样地优势高
    # 首先做判断：样地内每一株树都有树高时才计算
    data.Alive.h2 <- subset(data.h, data.h$d >= Dmin & !data.h$TreeType %in% c(13, 14, 15, 50))
    d_length <- tapply(data.Alive.h2$d, INDEX = data.Alive.h2$Plot, FUN = length)
    d_length <- data.frame("Plot" = rownames(as.data.frame(d_length)), d_length)

    data.Alive.h3 <- subset(data.h, data.h$d >= Dmin & data.h$h > 1.3 & !data.h$TreeType %in% c(13, 14, 15, 50))
    h_length <- tapply(data.Alive.h3$h, INDEX = data.Alive.h3$Plot, FUN = length)
    h_length <- data.frame("Plot" = rownames(as.data.frame(h_length)), h_length)

    temp_length <- left_join(d_length, h_length, by = "Plot")
    # 筛选出d和h具有相同行数的样地，计算优势高
    Plot_h <- temp_length %>%
      mutate(cc = d_length - h_length ) %>%
      filter(cc == 0) %>%
      select(Plot)

    if (nrow(Plot_h)!=0) {
      # 样地内最高的6株树的平均值， Ht
      # 样地内胸径最粗的6株树的平均值， Hd
      # 样地面积S每增加0.01ha，统计株数+1
      # 纯林不考虑树种
      # 混交林是否考虑树种--待商议
      Plot_h <- as.data.frame(Plot_h)
      colnames(Plot_h) <- "Plot"
      data.Alive.htStat <- left_join(as.data.frame(Plot_h), data.h, by = "Plot")

      h_TreeCounts <- round(S * 100, 0) # 根据样地面积确定的计算的株数

      # 样地内最高的 h_TreeCounts 株树的平均值
      # data.Alive.h1 = arrange(data.Alive.h1, Plot, h)
      data.Alive.ht <- data.Alive.htStat %>%
        group_by(Plot) %>%
        slice_max(n = h_TreeCounts, order_by = h, with_ties = F)
      dominantHt.Counts <- tapply(data.Alive.ht$h, INDEX = data.Alive.ht$Plot, FUN = length)
      dominantHt <- tapply(data.Alive.ht$h, INDEX = data.Alive.ht$Plot, FUN = mean)
      dominantHt[is.na(dominantHt)] <- 0
      dominantHt <- round(dominantHt, 1)
      dominantHt <- data.frame("Plot" = rownames(dominantHt), "Ht" = dominantHt, "Ht_Counts" = dominantHt.Counts)
      Alive.StandFactor3 <- left_join(Alive.StandFactor3, dominantHt, by = "Plot")

      # 样地内胸径最粗 h_TreeCounts 株树的平均值
      # data.Alive.h2 = arrange(data.Alive.h1, Plot, d)
      data.Alive.hd <- data.Alive.htStat %>%
        group_by(Plot) %>%
        arrange(desc(d), desc(h)) %>%
        slice_max(n = h_TreeCounts, order_by = d, with_ties = F)
      dominantHd.Counts <- tapply(data.Alive.hd$h, INDEX = data.Alive.hd$Plot, FUN = length)
      dominantHd <- tapply(data.Alive.hd$h, INDEX = data.Alive.hd$Plot, FUN = mean)
      dominantHd[is.na(dominantHd)] <- 0
      dominantHd <- round(dominantHd, 1)
      dominantHd <- data.frame("Plot" = rownames(dominantHd), "Hd" = dominantHd, "Hd_Counts" = dominantHd.Counts)
      Alive.StandFactor3 <- left_join(Alive.StandFactor3, dominantHd, by = "Plot")
    } else {
      cat("The dominant height (Ht and Hd) will not be calculated, because not every tree has tree height.\n")
    }
  }

  # 统计各样地的的蓄积、生物量、碳储量结果结果与上述林分因子结果合并
  ##
  # 统计蓄积
  ##
  if (!is.null(V)) {

    data <- data.frame("Plot" = Plot, "SP" = SP, "TreeType" = TT, "d" = D, "ba" = ba, "V" = V)
    data.Alive.V <- subset(data, data$d >= Dmin & !data$TreeType %in% c(13, 14, 15, 50))
    Alive.StandV <- tapply(data.Alive.V$V / S, INDEX = data.Alive.V$Plot, FUN = sum)
    Alive.StandV <- data.frame("Plot" = rownames(Alive.StandV), "V" = Alive.StandV)
    Alive.StandFactor3 <- left_join(Alive.StandFactor3, Alive.StandV, by = "Plot")

    if (length(unique(SP)) != 1) {
      # 蓄积比例----单独作为一个文件比较好，方便数据合并
      V.sum <- tapply(data.Alive.V$V, INDEX = list(data.Alive.V$Plot, data.Alive.V$SP), FUN = sum)
      V.sum[is.na(V.sum)] <- 0
      V.zhanbi <- round(V.sum / rowSums(V.sum) * 100, 0)
      V.zhanbi[is.na(V.zhanbi)] <- 0
      V.zhanbi[sapply(V.zhanbi, is.infinite)] <- 0
      # 归一化
      for (i in 1:nrow(V.zhanbi)) {
        zuchenghe <- sum(V.zhanbi[i, ])
        if (zuchenghe != 100) {
          if (zuchenghe < 100) {
            V.zhanbi[i, which.max(V.zhanbi[i, ])] <- V.zhanbi[i, which.max(V.zhanbi[i, ])] + (100 - zuchenghe)
          } else {
            V.zhanbi[i, which.max(V.zhanbi[i, ])] <- V.zhanbi[i, which.max(V.zhanbi[i, ])] - (zuchenghe - 100)
          }
        }
      }
      V.names <- colnames(V.zhanbi)
      V.names <- paste0("V.p_", V.names)
      colnames(V.zhanbi) <- V.names
      V.zhanbi <- data.frame("Plot" = rownames(as.data.frame(V.sum)), V.zhanbi)
      # 蓄积树种组成式
      V.Composition <- CompositionStat(V.zhanbi, start = 5)
      V.zhanbi <- left_join(V.zhanbi, V.Composition, by = "Plot")
    }
  }
  ##
  # 统计生物量
  ##
  if (!is.null(Biomass)) {
    Biomass <- as.data.frame(apply(Biomass, MARGIN = 2, FUN = as.numeric))
    data <- data.frame("Plot" = Plot, "TreeType" = TT, "d" = D, "ba" = ba, Biomass)
    data.Alive.Bio <- subset(data, data$d >= Dmin & !data$TreeType %in% c(13, 14, 15, 50))
    if (ncol(as.data.frame(Biomass)) != 1) {
      Alive.StandBio_1 <- tapply(data.Alive.Bio[, 4 + 1] / 1000 / S, INDEX = data.Alive.Bio$Plot, FUN = sum)
      Alive.StandBio_1 <- data.frame("Plot" = rownames(Alive.StandBio_1), "Bio_1" = Alive.StandBio_1)
      for (i in 2:ncol(Biomass)) {
        Alive.StandBio_i <- tapply(data.Alive.Bio[, 4 + i] / 1000 / S, INDEX = data.Alive.Bio$Plot, FUN = sum)
        Alive.StandBio_i <- data.frame(Alive.StandBio_i)
        colnames(Alive.StandBio_i) <- paste0("Bio_", i)
        Alive.StandBio_1 <- cbind(Alive.StandBio_1, Alive.StandBio_i)
      }
      Alive.StandFactor3 <- left_join(Alive.StandFactor3, Alive.StandBio_1, by = "Plot")
    } else {
      Alive.StandBio_1 <- tapply(data.Alive.Bio[, 4 + 1] / 1000 / S, INDEX = data.Alive.Bio$Plot, FUN = sum)
      Alive.StandBio_1 <- data.frame("Plot" = rownames(Alive.StandBio_1), "Bio" = Alive.StandBio_1)
      Alive.StandFactor3 <- left_join(Alive.StandFactor3, Alive.StandBio_1, by = "Plot")
    }
  }
  ##
  # 统计碳储量
  ##
  if (!is.null(Carbon)) {
    Carbon <- as.data.frame(apply(Carbon, MARGIN = 2, FUN = as.numeric))
    data <- data.frame("Plot" = Plot, "TreeType" = TT, "d" = D, "ba" = ba, Carbon)
    data.Alive.Car <- subset(data, data$d >= Dmin & !data$TreeType %in% c(13, 14, 15, 50))
    if (ncol(as.data.frame(Carbon)) != 1) {
      Alive.StandCar_1 <- tapply(data.Alive.Car[, 4 + 1] / 1000 / S, INDEX = data.Alive.Car$Plot, FUN = sum)
      Alive.StandCar_1 <- data.frame("Plot" = rownames(Alive.StandCar_1), "Car_1" = Alive.StandCar_1)
      for (i in 2:ncol(Carbon)) {
        Alive.StandCar_i <- tapply(data.Alive.Car[, 4 + i] / 1000 / S, INDEX = data.Alive.Car$Plot, FUN = sum)
        Alive.StandCar_i <- data.frame(Alive.StandCar_i)
        colnames(Alive.StandCar_i) <- paste0("Car_", i)
        Alive.StandCar_1 <- cbind(Alive.StandCar_1, Alive.StandCar_i)
      }
      Alive.StandFactor3 <- left_join(Alive.StandFactor3, Alive.StandCar_1, by = "Plot")
    } else {
      Alive.StandCar_1 <- tapply(data.Alive.Car[, 4 + 1] / 1000 / S, INDEX = data.Alive.Car$Plot, FUN = sum)
      Alive.StandCar_1 <- data.frame("Plot" = rownames(Alive.StandCar_1), "Car" = Alive.StandCar_1)
      Alive.StandFactor3 <- left_join(Alive.StandFactor3, Alive.StandCar_1, by = "Plot")
    }
  }

  # Cut_木和Dead_木的计算
  ##
  # Cut_木 条件：检尺类型%in%c(14,50))。
  # 50为非正常样木，2021-9 吉林5期数据单木逻辑检查时前期存在后期无记录的样木标记为非正常采伐木50
  ##
  if (sum(TreeType %in% c(14, 50)) != 0) {
    data <- data.frame("Plot" = Plot, "TreeType" = TT, "d" = D, "ba" = ba)
    data.Cut <- subset(data, data$d >= Dmin & data$TreeType %in% c(14, 50))
    if (nrow(data.Cut) != 0) {
      Cut.StandFactor <- tapply(data.Cut$d, INDEX = data.Cut$Plot, FUN = tempFUN)
      Cut.StandFactor2 <- data.frame(matrix(unlist(Cut.StandFactor, use.names = T), ncol = 5, byrow = T), stringsAsFactors = F)
      colnames(Cut.StandFactor2) <- c("Cut_Dg", "Cut_Ba", "Cut_N", "Cut_SDI", "Cut_aSDI")
      Cut.StandFactor3 <- data.frame("Plot" = rownames(Cut.StandFactor), Cut.StandFactor2)
      # 拼接
      Alive.StandFactor3 <- left_join(Alive.StandFactor3, Cut.StandFactor3, by = "Plot")
      # 统计Cut_木的蓄积  # 2021-10-16， 补充计算采伐木的蓄积
      if (!is.null(V)) {
        data <- data.frame("Plot" = Plot, "TreeType" = TT, "d" = D, "ba" = ba, "V" = V)
        data.Cut.V <- subset(data, data$d >= Dmin & data$TreeType %in% c(14, 50))
        Cut.StandV <- tapply(data.Cut.V$V / S, INDEX = data.Cut.V$Plot, FUN = sum)
        Cut.StandV <- data.frame("Plot" = rownames(Cut.StandV), "Cut_V" = Cut.StandV)
        # 拼接
        Alive.StandFactor3 <- left_join(Alive.StandFactor3, Cut.StandV, by = "Plot")
      }
      # 统计Cut_木的生物量 # 2021-10-17， 补充计算采伐木的生物量
      if (!is.null(Biomass)) {
        Biomass <- as.data.frame(apply(Biomass, MARGIN = 2, FUN = as.numeric))
        data <- data.frame("Plot" = Plot, "TreeType" = TT, "d" = D, "ba" = ba, Biomass)
        data.Cut.Bio <- subset(data, data$d >= Dmin & data$TreeType %in% c(14, 50))
        if (nrow(data.Cut.Bio) !=0){
          if (ncol(as.data.frame(Biomass)) != 1) {
            Cut.StandBio_1 <- tapply(data.Cut.Bio[, 4 + 1] / 1000 / S, INDEX = data.Cut.Bio$Plot, FUN = sum)
            Cut.StandBio_1 <- data.frame("Plot" = rownames(Cut.StandBio_1), "Cut_Bio_1" = Cut.StandBio_1)
            for (i in 2:ncol(Biomass)) {
              Cut.StandBio_i <- tapply(data.Cut.Bio[, 4 + i] / 1000 / S, INDEX = data.Cut.Bio$Plot, FUN = sum)
              Cut.StandBio_i <- data.frame(Cut.StandBio_i)
              colnames(Cut.StandBio_i) <- paste0("Cut_Bio_", i)
              Cut.StandBio_1 <- cbind(Cut.StandBio_1, Cut.StandBio_i)
            }
            Alive.StandFactor3 <- left_join(Alive.StandFactor3, Cut.StandBio_1, by = "Plot")
          } else {
            Cut.StandBio_1 <- tapply(data.Cut.Bio[, 4 + 1] / 1000 / S, INDEX = data.Cut.Bio$Plot, FUN = sum)
            Cut.StandBio_1 <- data.frame("Plot" = rownames(Cut.StandBio_1), "Cut_Bio" = Cut.StandBio_1)
            Alive.StandFactor3 <- left_join(Alive.StandFactor3, Cut.StandBio_1, by = "Plot")
          }
         }
        }


      # 统计Cut_木的碳储量 # 2021-10-17， 补充计算采伐木的碳储量
      if (!is.null(Carbon)) {
        Carbon <- as.data.frame(apply(Carbon, MARGIN = 2, FUN = as.numeric))
        data <- data.frame("Plot" = Plot, "TreeType" = TT, "d" = D, "ba" = ba, Carbon)
        data.Cut.Car <- subset(data, data$d >= Dmin & data$TreeType %in% c(14, 50))
        if (nrow(data.Cut.Car) !=0){
          if (ncol(as.data.frame(Carbon)) != 1) {
            Cut.StandCar_1 <- tapply(data.Cut.Car[, 4 + 1] / 1000 / S, INDEX = data.Cut.Car$Plot, FUN = sum)
            Cut.StandCar_1 <- data.frame("Plot" = rownames(Cut.StandCar_1), "Cut_Car_1" = Cut.StandCar_1)
            for (i in 2:ncol(Carbon)) {
              Cut.StandCar_i <- tapply(data.Cut.Car[, 4 + i] / 1000 / S, INDEX = data.Cut.Car$Plot, FUN = sum)
              Cut.StandCar_i <- data.frame(Cut.StandCar_i)
              colnames(Cut.StandCar_i) <- paste0("Cut_Car_", i)
              Cut.StandCar_1 <- cbind(Cut.StandCar_1, Cut.StandCar_i)
            }
            Alive.StandFactor3 <- left_join(Alive.StandFactor3, Cut.StandCar_1, by = "Plot")
          } else {
            Cut.StandCar_1 <- tapply(data.Cut.Car[, 4 + 1] / 1000 / S, INDEX = data.Cut.Car$Plot, FUN = sum)
            Cut.StandCar_1 <- data.frame("Plot" = rownames(Cut.StandCar_1), "Cut_Car" = Cut.StandCar_1)
            Alive.StandFactor3 <- left_join(Alive.StandFactor3, Cut.StandCar_1, by = "Plot")
          }
        }
      }
    }
  }
  ##
  # Dead_木 条件：检尺类型%in%c(13,15)
  ##
  if (sum(TreeType %in% c(13, 15)) != 0) {
    data <- data.frame("Plot" = Plot, "TreeType" = TT, "d" = D, "ba" = ba)
    data.Dead <- subset(data, data$d >= Dmin & data$TreeType %in% c(13, 15))
    if (nrow(data.Dead) != 0) {
      Dead.StandFactor <- tapply(data.Dead$d, INDEX = data.Dead$Plot, FUN = tempFUN)
      Dead.StandFactor2 <- data.frame(matrix(unlist(Dead.StandFactor, use.names = T), ncol = 5, byrow = T), stringsAsFactors = F)
      colnames(Dead.StandFactor2) <- c("Dead_Dg", "Dead_Ba", "Dead_N", "Dead_SDI", "Dead_aSDI")
      Dead.StandFactor3 <- data.frame("Plot" = rownames(Dead.StandFactor), Dead.StandFactor2)
      # 拼接
      Alive.StandFactor3 <- left_join(Alive.StandFactor3, Dead.StandFactor3, by = "Plot")
      # 统计Dead_木的蓄积   # 2021-10-16， 补充计算枯死木的蓄积
      if (!is.null(V)) {
        V <- as.numeric(V)
        data <- data.frame("Plot" = Plot, "TreeType" = TT, "d" = D, "ba" = ba, "V" = V)
        data.Dead.V <- subset(data, data$d >= Dmin & data$TreeType %in% c(13, 15))
        Dead.StandV <- tapply(data.Dead.V$V / S, INDEX = data.Dead.V$Plot, FUN = sum)
        Dead.StandV <- data.frame("Plot" = rownames(Dead.StandV), "Dead_V" = Dead.StandV)
        # 拼接
        Alive.StandFactor3 <- left_join(Alive.StandFactor3, Dead.StandV, by = "Plot")
      }
      # 统计Dead_木的生物量 # 2021-10-17， 补充计算枯死木的生物量
      if (!is.null(Biomass)) {
        Biomass <- as.data.frame(apply(Biomass, MARGIN = 2, FUN = as.numeric))
        data <- data.frame("Plot" = Plot, "TreeType" = TT, "d" = D, "ba" = ba, Biomass)
        data.Dead.Bio <- subset(data, data$d >= Dmin & data$TreeType %in% c(13, 15))
        if (nrow(data.Dead.Bio)!=0){
          if (ncol(as.data.frame(Biomass)) != 1) {
            Dead.StandBio_1 <- tapply(data.Dead.Bio[, 4 + 1] / 1000 / S, INDEX = data.Dead.Bio$Plot, FUN = sum)
            Dead.StandBio_1 <- data.frame("Plot" = rownames(Dead.StandBio_1), "Dead_Bio_1" = Dead.StandBio_1)
            for (i in 2:ncol(Biomass)) {
              Dead.StandBio_i <- tapply(data.Dead.Bio[, 4 + i] / 1000 / S, INDEX = data.Dead.Bio$Plot, FUN = sum)
              Dead.StandBio_i <- data.frame(Dead.StandBio_i)
              colnames(Dead.StandBio_i) <- paste0("Dead_Bio_", i)
              Dead.StandBio_1 <- cbind(Dead.StandBio_1, Dead.StandBio_i)
            }
            Alive.StandFactor3 <- left_join(Alive.StandFactor3, Dead.StandBio_1, by = "Plot")
          } else {
            Dead.StandBio_1 <- tapply(data.Dead.Bio[, 4 + 1] / 1000 / S, INDEX = data.Dead.Bio$Plot, FUN = sum)
            Dead.StandBio_1 <- data.frame("Plot" = rownames(Dead.StandBio_1), "Dead_Bio" = Dead.StandBio_1)
            Alive.StandFactor3 <- left_join(Alive.StandFactor3, Dead.StandBio_1, by = "Plot")
          }
        }
      }
      # 统计Dead_木的碳储量 # 2021-10-17， 补充计算枯死木的碳储量
      if (!is.null(Carbon)) {
        Carbon <- as.data.frame(apply(Carbon, MARGIN = 2, FUN = as.numeric))
        data <- data.frame("Plot" = Plot, "TreeType" = TT, "d" = D, "ba" = ba, Carbon)
        data.Dead.Car <- subset(data, data$d >= Dmin & data$TreeType %in% c(13, 15))
        if (nrow(data.Dead.Car)!=0){
          if (ncol(as.data.frame(Carbon)) != 1) {
            Dead.StandCar_1 <- tapply(data.Dead.Car[, 4 + 1] / 1000 / S, INDEX = data.Dead.Car$Plot, FUN = sum)
            Dead.StandCar_1 <- data.frame("Plot" = rownames(Dead.StandCar_1), "Dead_Car_1" = Dead.StandCar_1)
            for (i in 2:ncol(Carbon)) {
              Dead.StandCar_i <- tapply(data.Dead.Car[, 4 + i] / 1000 / S, INDEX = data.Dead.Car$Plot, FUN = sum)
              Dead.StandCar_i <- data.frame(Dead.StandCar_i)
              colnames(Dead.StandCar_i) <- paste0("Dead_Car_", i)
              Dead.StandCar_1 <- cbind(Dead.StandCar_1, Dead.StandCar_i)
            }
            Alive.StandFactor3 <- left_join(Alive.StandFactor3, Dead.StandCar_1, by = "Plot")
          } else {
            Dead.StandCar_1 <- tapply(data.Dead.Car[, 4 + 1] / 1000 / S, INDEX = data.Dead.Car$Plot, FUN = sum)
            Dead.StandCar_1 <- data.frame("Plot" = rownames(Dead.StandCar_1), "Dead_Car" = Dead.StandCar_1)
            Alive.StandFactor3 <- left_join(Alive.StandFactor3, Dead.StandCar_1, by = "Plot")
          }
        }
      }
    }
  }

  # 其他类型木，如健康状态、上下层、目标树和非目标树等不同水平的林分因子
  if (!is.null(TreeClass) & length(unique(TreeClass)) == 1) {
    cat("The 'TreeClass' was only one level, please chack it!\n")
  } else{
    TreeClass <- as.character(TreeClass)
    if(sum(is.na(TreeClass))!=0){
      cat(paste0(Na.f(TreeClass), "% Missing value (or NA) in 'TreeClass'.\n"))
    }
  }
  if (!is.null(TreeClass) & length(unique(TreeClass)) > 1) {
    data <- data.frame("Plot" = Plot, "TreeClass" = paste0("TreeClass", TreeClass), "TreeType" = TT, "d" = D)
    data.Alive <- subset(data, (data$d >= Dmin) & (!data$TreeType %in% c(13, 14, 15, 50)))
    # Dg
    Alive.TreeClass.Dg <- tapply(data.Alive$d, INDEX = list(data.Alive$Plot, data.Alive$TreeClass), FUN = Dg)
    Alive.TreeClass.Dg[is.na(Alive.TreeClass.Dg)] <- 0
    Alive.TreeClass.Dg.names <- colnames(Alive.TreeClass.Dg)
    Alive.TreeClass.Dg.names <- paste0(Alive.TreeClass.Dg.names, "_Dg")
    colnames(Alive.TreeClass.Dg) <- Alive.TreeClass.Dg.names
    # Ba
    Alive.TreeClass.Ba <- tapply(data.Alive$d, INDEX = list(data.Alive$Plot, data.Alive$TreeClass), FUN = Ba)
    Alive.TreeClass.Ba[is.na(Alive.TreeClass.Ba)] <- 0
    Alive.TreeClass.Ba.names <- colnames(Alive.TreeClass.Ba)
    Alive.TreeClass.Ba.names <- paste0(Alive.TreeClass.Ba.names, "_Ba")
    colnames(Alive.TreeClass.Ba) <- Alive.TreeClass.Ba.names
    # N
    Alive.TreeClass.N <- tapply(data.Alive$d, INDEX = list(data.Alive$Plot, data.Alive$TreeClass), FUN = N)
    Alive.TreeClass.N[is.na(Alive.TreeClass.N)] <- 0
    Alive.TreeClass.N.names <- colnames(Alive.TreeClass.N)
    Alive.TreeClass.N.names <- paste0(Alive.TreeClass.N.names, "_N")
    colnames(Alive.TreeClass.N) <- Alive.TreeClass.N.names
    # SDI
    Alive.TreeClass.SDI <- tapply(data.Alive$d, INDEX = list(data.Alive$Plot, data.Alive$TreeClass), FUN = SDI)
    Alive.TreeClass.SDI[is.na(Alive.TreeClass.SDI)] <- 0
    Alive.TreeClass.SDI.names <- colnames(Alive.TreeClass.SDI)
    Alive.TreeClass.SDI.names <- paste0(Alive.TreeClass.SDI.names, "_SDI")
    colnames(Alive.TreeClass.SDI) <- Alive.TreeClass.SDI.names
    # aSDI
    Alive.TreeClass.aSDI <- tapply(data.Alive$d, INDEX = list(data.Alive$Plot, data.Alive$TreeClass), FUN = aSDI)
    Alive.TreeClass.aSDI[is.na(Alive.TreeClass.aSDI)] <- 0
    Alive.TreeClass.aSDI.names <- colnames(Alive.TreeClass.aSDI)
    Alive.TreeClass.aSDI.names <- paste0(Alive.TreeClass.aSDI.names, "_aSDI")
    colnames(Alive.TreeClass.aSDI) <- Alive.TreeClass.aSDI.names

    # 分林木类型的基本林分因子
    Alive.TreeClass <- data.frame("Plot" = rownames(Alive.TreeClass.Dg), Alive.TreeClass.Dg, Alive.TreeClass.Ba, Alive.TreeClass.N, Alive.TreeClass.SDI, Alive.TreeClass.aSDI, check.names = F)
    # 拼接
    Alive.StandFactor3 <- left_join(Alive.StandFactor3, Alive.TreeClass, by = "Plot")

    # 分林木类型的蓄积  # 2021-10-16， 补充计算其他类型木的蓄积
    if (!is.null(V)) {
      V <- as.numeric(V)
      data <- data.frame("Plot" = Plot, "TreeClass" = paste0("TreeClass", TreeClass), "TreeType" = TT, "d" = D, "V" = V)
      data.Alive.TreeClass.V <- subset(data, (data$d >= Dmin) & (!data$TreeType %in% c(13, 14, 15, 50)))
      Alive.TreeClass.V <- tapply(data.Alive.TreeClass.V$V / S, INDEX = list(data.Alive.TreeClass.V$Plot, data.Alive.TreeClass.V$TreeClass), FUN = sum)
      Alive.TreeClass.V[is.na(Alive.TreeClass.V)] <- 0
      Alive.TreeClass.V.names <- colnames(Alive.TreeClass.V)
      Alive.TreeClass.V.names <- paste0(Alive.TreeClass.V.names, "_V")
      colnames(Alive.TreeClass.V) <- Alive.TreeClass.V.names
      Alive.TreeClass.V <- data.frame("Plot" = rownames(Alive.TreeClass.V), Alive.TreeClass.V, check.name = F)
      # 拼接
      Alive.StandFactor3 <- left_join(Alive.StandFactor3, Alive.TreeClass.V, by = "Plot")
    }
    # 分林木类型的生物量  # 2021-10-16， 补充计算其他类型木的生物量
    if (!is.null(Biomass)) {
      Biomass <- as.data.frame(apply(Biomass, MARGIN = 2, FUN = as.numeric))
      data <- data.frame("Plot" = Plot, "TreeClass" = paste0("TreeClass", TreeClass), "TreeType" = TT, "d" = D, Biomass)
      data.TreeClass.Bio <- subset(data, (data$d >= Dmin) & (!data$TreeType %in% c(13, 14, 15, 50)))
      if (nrow(data.TreeClass.Bio)!=0){
        if (ncol(as.data.frame(Biomass)) != 1) {
          TreeClass.StandBio_1 <- tapply(data.TreeClass.Bio[, 4 + 1] / 1000 / S, INDEX = list(data.TreeClass.Bio$Plot, data.TreeClass.Bio$TreeClass), FUN = sum)
          TreeClass.StandBio_1[is.na(TreeClass.StandBio_1)] <- 0
          TreeClass.StandBio_1.names <- colnames(TreeClass.StandBio_1)
          TreeClass.StandBio_1.names <- paste0(TreeClass.StandBio_1.names, "_Bio_1")
          colnames(TreeClass.StandBio_1) <- TreeClass.StandBio_1.names
          TreeClass.StandBio_1 <- data.frame("Plot" = rownames(TreeClass.StandBio_1), TreeClass.StandBio_1, check.name = F)
          for (i in 2:ncol(Biomass)) {
            TreeClass.StandBio_i <- tapply(data.TreeClass.Bio[, 4 + i] / 1000 / S, INDEX = list(data.TreeClass.Bio$Plot, data.TreeClass.Bio$TreeClass), FUN = sum)
            TreeClass.StandBio_i[is.na(TreeClass.StandBio_i)] <- 0
            TreeClass.StandBio_i.names <- colnames(TreeClass.StandBio_i)
            colnames(TreeClass.StandBio_i) <- paste0(TreeClass.StandBio_i.names, "_Bio_", i)
            TreeClass.StandBio_1 <- cbind(TreeClass.StandBio_1, TreeClass.StandBio_i)
          }
          Alive.StandFactor3 <- left_join(Alive.StandFactor3, TreeClass.StandBio_1, by = "Plot")
        } else {
          TreeClass.StandBio_1 <- tapply(data.TreeClass.Bio[, 4 + 1] / 1000 / S, INDEX = list(data.TreeClass.Bio$Plot, data.TreeClass.Bio$TreeClass), FUN = sum)
          TreeClass.StandBio_1[is.na(TreeClass.StandBio_1)] <- 0
          TreeClass.StandBio_1.names <- colnames(TreeClass.StandBio_1)
          TreeClass.StandBio_1.names <- paste0(TreeClass.StandBio_1.names, "_Bio")
          colnames(TreeClass.StandBio_1) <- TreeClass.StandBio_1.names
          TreeClass.StandBio_1 <- data.frame("Plot" = rownames(TreeClass.StandBio_1), TreeClass.StandBio_1, check.name = F)
          Alive.StandFactor3 <- left_join(Alive.StandFactor3, TreeClass.StandBio_1, by = "Plot")
        }
      }
    }
    # 分林木类型的碳储量  # 2021-10-16， 补充计算其他类型木的碳储量
    if (!is.null(Carbon)) {
      Carbon <- as.data.frame(apply(Carbon, MARGIN = 2, FUN = as.numeric))
      data <- data.frame("Plot" = Plot, "TreeClass" = paste0("TreeClass", TreeClass), "TreeType" = TT, "d" = D, Carbon)
      data.TreeClass.Car <- subset(data, (data$d >= Dmin) & (!data$TreeType %in% c(13, 14, 15, 50)))
      if (nrow(data.TreeClass.Car)!=0){
        if (ncol(as.data.frame(Carbon)) != 1) {
          TreeClass.StandCar_1 <- tapply(data.TreeClass.Car[, 4 + 1] / 1000 / S, INDEX = list(data.TreeClass.Car$Plot, data.TreeClass.Car$TreeClass), FUN = sum)
          TreeClass.StandCar_1[is.na(TreeClass.StandCar_1)] <- 0
          TreeClass.StandCar_1.names <- colnames(TreeClass.StandCar_1)
          TreeClass.StandCar_1.names <- paste0(TreeClass.StandCar_1.names, "_Car_1")
          colnames(TreeClass.StandCar_1) <- TreeClass.StandCar_1.names
          TreeClass.StandCar_1 <- data.frame("Plot" = rownames(TreeClass.StandCar_1), TreeClass.StandCar_1)
          for (i in 2:ncol(Carbon)) {
            TreeClass.StandCar_i <- tapply(data.TreeClass.Car[, 4 + i] / 1000 / S, INDEX = list(data.TreeClass.Car$Plot, data.TreeClass.Car$TreeClass), FUN = sum)
            TreeClass.StandCar_i[is.na(TreeClass.StandCar_i)] <- 0
            TreeClass.StandCar_i.names <- colnames(TreeClass.StandCar_i)
            colnames(TreeClass.StandCar_i) <- paste0(TreeClass.StandCar_i.names, "_Car_", i)
            TreeClass.StandCar_1 <- cbind(TreeClass.StandCar_1, TreeClass.StandCar_i)
          }
          Alive.StandFactor3 <- left_join(Alive.StandFactor3, TreeClass.StandCar_1, by = "Plot")
        } else {
          TreeClass.StandCar_1 <- tapply(data.TreeClass.Car[, 4 + 1] / 1000 / S, INDEX = list(data.TreeClass.Car$Plot, data.TreeClass.Car$TreeClass), FUN = sum)
          TreeClass.StandCar_1[is.na(TreeClass.StandCar_1)] <- 0
          TreeClass.StandCar_1.names <- colnames(TreeClass.StandCar_1)
          TreeClass.StandCar_1.names <- paste0(TreeClass.StandCar_1.names, "_Car")
          colnames(TreeClass.StandCar_1) <- TreeClass.StandCar_1.names
          TreeClass.StandCar_1 <- data.frame("Plot" = rownames(TreeClass.StandCar_1), TreeClass.StandCar_1)
          Alive.StandFactor3 <- left_join(Alive.StandFactor3, TreeClass.StandCar_1, by = "Plot")
        }
      }
    }
  }

  Alive.StandFactor3[is.na(Alive.StandFactor3)] <- 0

  rownames(Alive.StandFactor3) <- NULL

  # 控制输出
  if (!is.null(V) & length(unique(SP)) != 1) {
    rownames(ba.zhanbi) <- NULL
    rownames(n.zhanbi) <- NULL
    rownames(V.zhanbi) <- NULL
    rownames(Alive.SP) <- NULL
    list.result <- list("Stand Factors" = Alive.StandFactor3, "Ba Composition" = ba.zhanbi, "N Composition" = n.zhanbi, "V Composition" = V.zhanbi, "Species Stand Factors" = Alive.SP)
  } else if (length(unique(SP)) != 1) {
    rownames(ba.zhanbi) <- NULL
    rownames(n.zhanbi) <- NULL
    rownames(Alive.SP) <- NULL
    list.result <- list("Stand Factors" = Alive.StandFactor3, "Ba Composition" = ba.zhanbi, "N Composition" = n.zhanbi, "Species Stand Factors" = Alive.SP)
  } else {
    list.result <- list("Stand Factors" = Alive.StandFactor3)
  }

  return(list.result)
}
