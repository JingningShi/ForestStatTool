#' @title 重要值计算

#' @author Xiao He: hexiaonuist@163.com
#' @author Chaofan Zhou: cfzhou2021@163.com
#' @author Guangshuang Duan: oliverdgs@163.com
#' @author Jingning Shi: shijingning@hebau.edu.cn
#' @author Xiangdong Lei: xdlei@ifrit.ac.cn

#' @usage IV (D, Plot, SubPlot, SP, TreeType = NULL, Dmin = 5)

#' @description 此函数实现了输入单木测树因子统计输出【相对频度 Fr】、【相对丰度 Dr】、【对优势度 Pr】和【重要值IV】

#' @param D （必须）数值向量。样地林木胸径。单位：cm。
#' @param Dmin 数值。起测胸径，默认为5cm。
#' @param Plot （必须）向量，可以是数值型，也可以是字符型。样地号。
#' @param SubPlot （必须）向量，可以是数值型，也可以是字符型。子样地号。
#' @param SP （必须）向量，可以是数值型，也可以是字符型。树种。
#' @param TreeType （可选）数值向量。检尺类型。12和14表示枯死木，13和50表示采伐木，其余值均认为是活立木，可输入全部检尺类型，程序将筛选出活立木（检尺类型不等于12、13、14和50）用于计算重要值。默认为空，即输入的数据均为活立木。

#' @return 1）输出结果为一个数据框组成
#' @return 2）输出变量名：
#' @return Plot：样地号
#' @return Fr_XXX：XXX树种的相对频度
#' @return Dr_XXX：XXX树种的相对丰度
#' @return Pr_XXX：XXX树种的相对优势度
#' @return IV_XXX：XXX树种的重要值

#' @export IV
#' @name IV

#' @examples ## 加载内置数据
#' @examples data(ForestStatTool)
#' @examples #bigplot是一个100m×100m的矩形样地的林木信息，首先利用Plot_Divide函数将bigplot划分为后的25个20m×20m的小样地的合并数据框。
#' @examples #除原始变量外，还包含的新列有：分割后的新样地名subplot，小样地的新X坐标subX和小样地的新Y坐标subY
#' @examples subplot1 <- Plot_Divide(Data = bigplot, Plot = bigplot$Plot, X = bigplot$X, Y = bigplot$Y,
#' @examples                         Num_xy = c(5,5), Range_xy = c(100,100))
#' @examples #计算bigplot重要值
#' @examples a1 = IV(Plot = subplot1$Plot, SubPlot = subplot1$subplot, SP = subplot1$SP, D = subplot1$D)
#' @examples head(a1)


IV <- function(D, Plot, SubPlot, SP, TreeType = NULL, Dmin = 5) {
  # 统计NA的比例，何潇-2022-12-8
  Na.f=function(x){ round(sum(is.na(x))/length(x)*100, 2) }
  if(is.null(Plot)){
    stop("'Plot' don't input.")
  }else if(sum(is.na(Plot))!=0){
    print(paste0(Na.f(Plot), "% Missing value (or NA) in 'Plot'."))
  }
  if(is.null(SubPlot)){
    stop("'SubPlot' don't input.")
  }else if(sum(is.na(SubPlot))!=0){
    print(paste0(Na.f(SubPlot), "% Missing value (or NA) in 'SubPlot'."))
  }
  n1 = nlevels(factor(Plot))
  n2 = nlevels(factor(paste0(Plot,"_",SubPlot)))
  if(n1>n2){
    stop("The number of 'SubPlot' is lower than 'Plot'.")
  }
  if(is.null(D)|!is.numeric(D)){
    stop("'D' don't input or 'D' must be numeric.")
  }else if(sum(is.na(D))!=0){
    print(paste0(Na.f(D), "% Missing value (or NA) in 'D'."))
  }
  if(is.null(SP)){
    stop("'SP' don't input.")
  }else if(sum(is.na(SP))!=0){
    print(paste0(Na.f(SP), "% Missing value (or NA) in 'SP'."))
  }

  Plot <- as.character(Plot)
  SubPlot <- as.character(SubPlot)
  ba <- pi * D^2 / 40000
  if (!is.null(TreeType)) {
    data <- data.frame("Plot" = Plot, "SubPlot" = SubPlot, "SP" = SP, "TreeType" = TreeType, "d" = D, "ba" = ba)
    data.Alive <- subset(data, data$d >= Dmin & !data$TreeType %in% c(13, 14, 15, 50))
  } else {
    data <- data.frame("Plot" = Plot, "SubPlot" = SubPlot, "SP" = SP, "d" = D, "ba" = ba)
    data.Alive <- subset(data, data$d >= Dmin)
  }

  # 相对频度 Fr
  number <- function(x) {
    length(unique(x))
  }
  number.sp <- tapply(data.Alive$SubPlot, INDEX = list(data.Alive$Plot, as.character(data.Alive$SP)), FUN = number)
  number.sp[is.na(number.sp)] <- 0
  Fr <- number.sp / rowSums(number.sp) * 100
  Fr <- as.data.frame(Fr)
  name.Fr <- colnames(Fr)
  name.Fr <- paste0("Fr_", name.Fr)
  colnames(Fr) <- name.Fr
  # 相对丰度 Dr
  counts.sp <- tapply(data.Alive$d, INDEX = list(data.Alive$Plot, as.character(data.Alive$SP)), FUN = length)
  counts.sp[is.na(counts.sp)] <- 0
  Dr <- counts.sp / rowSums(counts.sp) * 100
  Dr <- as.data.frame(Dr)
  name.Dr <- colnames(Dr)
  name.Dr <- paste0("Dr_", name.Dr)
  colnames(Dr) <- name.Dr
  # 相对优势度 Pr
  ba.sp <- tapply(data.Alive$ba, INDEX = list(data.Alive$Plot, as.character(data.Alive$SP)), FUN = sum)
  ba.sp[is.na(ba.sp)] <- 0
  Pr <- ba.sp / rowSums(ba.sp) * 100
  Pr <- as.data.frame(Pr)
  name.Pr <- colnames(Pr)
  name.Pr <- paste0("Pr_", name.Pr)
  colnames(Pr) <- name.Pr
  # 重要值IV
  IV <- (Fr + Dr + Pr) / 3
  IV <- as.data.frame(IV)
  name.IV <- colnames(IV)
  n <- nchar(name.IV)
  name.IV <- substring(name.IV, first = 4)
  name.IV <- paste0("IV_", name.IV)
  colnames(IV) <- name.IV
  IV <- data.frame("Plot" = rownames(as.data.frame(ba.sp)), Fr, Dr, Pr, IV)
  return(IV)
}
