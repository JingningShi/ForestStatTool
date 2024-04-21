#' @title  Tree species diversity indices. 树种多样性指数

#' @author Xiao He: hexiaonuist@163.com
#' @author Chaofan Zhou: cfzhou2021@163.com
#' @author Guangshuang Duan: oliverdgs@163.com
#' @author Jingning Shi: shijingning@hebau.edu.cn
#' @author Xiangdong Lei: xdlei@ifrit.ac.cn

#' @usage Div(Plot, D, SP, TreeType = NULL, Dmin = 5,  dClass = NULL, Index = 'Species')

#' @description 此函数实现了输入单木胸径和树种等调查变量，统计多样性指数。
#' @description 具体还包括物种多样性、（非空间）结构多样性等。

#' @param D （必须）数值向量。样地林木胸径。单位：cm。
#' @param Plot （必须）向量，可以是数值型，也可以是字符型。样地号。
#' @param SP （可选）向量，可以是数值型，也可以是字符型。树种。
#' @param TreeType （可选）数值向量。检尺类型。12和14表示枯死木，13和50表示采伐木，其余值均认为是活立木，可输入全部检尺类型，程序将自动计算活立木（检尺类型不等于12、13、14和50）的生物多样性因子。默认为空，即输入的数据均为活立木。
#' @param Dmin 数值。起测胸径，默认为5cm。
#' @param dClass 数值。径阶，2或4。其值为空时会同时输出2cm和4cm径阶划分所计算的结构多样性指数，默认为空。
#' @param Index 字符串。指标名称，取值为'Species'、'Size'、'SS'、'ALL_2'、'ALL_4'、'ALL'之一。'Species'，计算物种多样性/均匀度指数；'Size'，计算大小多样性/均匀度指数；'SS'，计算物种大小综合多样性/均匀度指数；搭配dClass参数可实现只计算2cm或只计算4cm径阶划分的大小多样性/均匀度指数、综合多样性/均匀度指数；'ALL_2'，计算物种多样性/均匀度指数、2cm径阶划分的大小多样性/均匀度指数以及2cm径阶划分的综合多样性/均匀度指数；'ALL_4'，计算物种多样性/均匀度指数、4cm径阶划分的大小多样性/均匀度指数以及4cm径阶划分的综合多样性/均匀度指数；'ALL'，计算所有多样性/均匀度指数。当Index取值为'ALL_2'、'ALL_4'、 'ALL'之一时，dClass参数不起作用。默认为'Species'。

#' @return 输出结果为一个数据框。
#' @return 2）输出变量名：
#' @return Plot：样地号
#' @return Richness：树种丰富度指数
#' @return N_GiniSimpson_Species_Div：Pi为株数占比时的GiniSimpson树种多样性指数
#' @return N_Simpson_Species_Div：Pi为株数占比时的Simpson树种多样性指数
#' @return N_Shannon_Species_Div：Pi为株数占比时的Shannon树种多样性指数
#' @return Margalef_Div：Margalef丰富度指数
#' @return Menhinick_Div：Menhinick丰富度指数
#' @return McIntosh_Div：McIntosh丰富度指数
#' @return BA_GiniSimpson_Species_Div：Pi为断面积占比时的GiniSimpson树种多样性指数
#' @return BA_Simpson_Species_Div：Pi为断面积占比时的Simpson树种多样性指数
#' @return BA_Shannon_Species_Div：Pi为断面积占比时的Shannon树种多样性指数
#' @return Chao1：Chao1指数
#' @return N_GiniSimpson_Size2_Div：Pi为株数占比时的GiniSimpson大小多样性指数（按2cm划分径阶）
#' @return N_Simpson_Size2_Div：Pi为株数占比时的Simpson大小多样性指数（按2cm划分径阶）
#' @return N_Shannon_Size2_Div：Pi为株数占比时的Shannon大小多样性指数（按2cm划分径阶）
#' @return BA_GiniSimpson_Size2_Div：Pi为断面积占比时的GiniSimpson大小多样性指数（按2cm划分径阶）
#' @return BA_Simpson_Size2_Div：Pi为断面积占比时的Simpson大小多样性指数（按2cm划分径阶）
#' @return BA_Shannon_Size2_Div：Pi为断面积占比时的Shannon大小多样性指数（按2cm划分径阶）
#' @return N_GiniSimpson_Size4_Div：Pi为株数占比时的GiniSimpson大小多样性指数（按4cm划分径阶）
#' @return N_Simpson_Size4_Div：Pi为株数占比时的Simpson大小多样性指数（按4cm划分径阶）
#' @return N_Shannon_Size4_Div：Pi为株数占比时的Shannon大小多样性指数（按4cm划分径阶）
#' @return BA_GiniSimpson_Size4_Div：Pi为断面积占比时的GiniSimpson大小多样性指数（按4cm划分径阶）
#' @return BA_Simpson_Size4_Div：Pi为断面积占比时的Simpson大小多样性指数（按4cm划分径阶）
#' @return BA_Shannon_Size4_Div：Pi为断面积占比时的Shannon大小多样性指数（按4cm划分径阶）
#' @return N_GiniSimpson_SS2_Div：Pi为株数占比时的GiniSimpson综合多样性指数（按2cm划分径阶）
#' @return N_Simpson_SS2_Div：Pi为株数占比时的Simpson综合多样性指数（按2cm划分径阶）
#' @return N_Shannon_SS2_Div：Pi为株数占比时的Shannon综合多样性指数（按2cm划分径阶）
#' @return BA_GiniSimpson_SS2_Div：Pi为断面积占比时的GiniSimpson综合多样性指数（按2cm划分径阶）
#' @return BA_Simpson_SS2_Div：Pi为断面积占比时的Simpson综合多样性指数（按2cm划分径阶）
#' @return BA_Shannon_SS2_Div：Pi为断面积占比时的Shannon综合多样性指数（按2cm划分径阶）
#' @return N_GiniSimpson_SS4_Div：Pi为株数占比时的GiniSimpson综合多样性指数（按4cm划分径阶）
#' @return N_Simpson_SS4_Div：Pi为株数占比时的Simpson综合多样性指数（按4cm划分径阶）
#' @return N_Shannon_SS4_Div：Pi为株数占比时的Shannon综合多样性指数（按4cm划分径阶）
#' @return BA_GiniSimpson_SS4_Div：Pi为断面积占比时的GiniSimpson综合多样性指数（按4cm划分径阶）
#' @return BA_Simpson_SS4_Div：Pi为断面积占比时的Simpson综合多样性指数（按4cm划分径阶）
#' @return BA_Shannon_SS4_Div：Pi为断面积占比时的Shannon综合多样性指数（按4cm划分径阶）

#' @export Div
#' @name Div

#' @import dplyr
#' @import stats

#' @examples ## 加载内置数据
#' @examples data(ForestStatTool)
#' @examples #bigplot是一个100m×100m的矩形样地的林木信息，首先利用Plot_Divide函数将bigplot划分为后的25个20m×20m的小样地的合并数据框。
#' @examples #除原始变量外，还包含的新列有：分割后的新样地名subplot，小样地的新X坐标subX和小样地的新Y坐标subY
#' @examples subplot1 <- Plot_Divide(Data = bigplot, Plot = bigplot$Plot, X = bigplot$X, Y = bigplot$Y,
#' @examples                       Num_xy = c(5,5), Range_xy = c(100,100))
#' @examples #计算不同小样地的多样性指标
#' @examples a1 = Div(Plot = subplot1$subplot, D = subplot1$D, SP = subplot1$SP, dClass = NULL, Index = 'ALL')
#' @examples head(a1)


Div = function(Plot, D, SP, TreeType = NULL, Dmin = 5,  dClass = NULL, Index = 'Species'){
  # 检测dplyr包是否存在，只有不存在时才会安装
  if(!require ('dplyr')){
    install.packages("dplyr")
    library(dplyr)
  }
  # 统计NA的比例，何潇-2022-12-8
  Na.f=function(x){ round(sum(is.na(x))/length(x)*100, 2) }
  # 数据组织
  data <- data.frame('Plot'=as.character(Plot), 'SP'=as.character(SP),  'd'=D)
  # 变量检查
  if(sum(is.na(data$Plot))!=0){
    cat(paste0(Na.f(data$Plot), "% Missing value (or NA) in 'Plot'.\n"))
  }
  if(!is.numeric(data$d)){
    stop("'D' must be numeric.")
  }else if(sum(is.na(data$d))!=0){
    cat(paste0(Na.f(data$d), "% Missing value (or NA) in 'D'.\n"))
  }
  if(!is.null(data$SP)){
    if(T %in% is.na(data$SP)){
      cat(paste0(Na.f(data$SP), "Missing value (or NA) in 'SP'.\n"))
    }
  }
  if(!is.null(TreeType)){
    data$TreeType = TreeType
    data.Alive = subset(data, data$d>=Dmin&!data$TreeType%in%c(13,14,15,50))
    }else{
    data.Alive = subset(data, data$d>=Dmin)
  }
  data.Alive$ba = pi*data.Alive$d^2/40000

  #删除NA值防止无法运算，何潇-2022-11-27
  data.Alive <- na.omit(data.Alive)
  if(nrow(data.Alive)==0){
    temp = apply(data.Alive, 2, function(x) sum(!is.na(x))==0)
    tempName = names(which(temp == TRUE))
    if(length(tempName)==1){
      stop(paste0("The column of", tempName, "is NA."))
    }else{
      stop(paste0("The columns of", tempName, "are NA."))
    }
  }

  ###################
  speciesdiv = function(comm, method = "full", tol = 1e-08){
    if (any(!method %in% c("richness", "GiniSimpson",
                           "Simpson", "Shannon", "Margalef", "Menhinick",
                           "McIntosh", "full")))
      stop("Your choice for method is not available")
    if ("full" %in% method)
      method <- c("richness", "GiniSimpson", "Simpson",
                  "Shannon", "Margalef", "Menhinick",
                  "McIntosh")
    if (any(comm < (-tol)))
      stop("Abundance entry in comm must be nonnegative")
    comm[comm < tol] <- 0
    if (all(rowSums(comm) < tol))
      stop("All communities are empty")
    if (any(rowSums(comm) < tol))
      warning("Empty communities were discarded")
    comm <- comm[rowSums(comm) > tol, , drop = FALSE]
    FUNshannon <- function(v) {
      if (length(v[v > 0]) == 1)
        return(0)
      else {
        v <- v[v > 0]
        return(as.vector(-sum(v/sum(v) * log(v/sum(v)))))
      }
    }
    RES <- matrix(0, nrow(comm), length(method))
    rownames(RES) <- rownames(comm)
    colnames(RES) <- method
    for (i in 1:length(method)) {
      if (method[i] == "richness")
        RES[, i] <- apply(comm, 1, function(x) length(x[x > 0]))
      if (method[i] == "GiniSimpson")
        RES[, i] <- apply(comm, 1, function(x) 1 - sum((x/sum(x))^2))
      if (method[i] == "Simpson")
        RES[, i] <- apply(comm, 1, function(x) 1/sum((x/sum(x))^2))
      if (method[i] == "Shannon")
        RES[, i] <- apply(comm, 1, FUNshannon)
      if (method[i] == "Margalef")
        RES[, i] <- apply(comm, 1, function(x) (length(x[x > 0]) - 1)/log(sum(x)))
      if (method[i] == "Menhinick")
        RES[, i] <- apply(comm, 1, function(x) length(x[x > 0])/sqrt(sum(x)))
      if (method[i] == "McIntosh")
        RES[, i] <- apply(comm, 1, function(x) (sum(x) - sqrt(sum(x^2)))/(sum(x) - sqrt(sum(x))))
    }
    return(RES)
  }
  #######################

  ####
  ## 物种 Species
  ## 包括物种丰富度
  ## 以及adiv包函数中speciesdiv()[6个指数]的计算结果，可基于株数n或断面积ba计算。以及Chao1指数。 # Chao1指数，何潇添加，2024-4-7
  ####
  # 基于株数n
  # 多样性指数
  counts.sp = tapply(data.Alive$d,INDEX=list(data.Alive$Plot, data.Alive$SP), FUN=length)
  counts.sp[is.na(counts.sp)]=0
  Ndiv = speciesdiv(counts.sp)
  Ndiv_out = data.frame("Plot"=rownames(as.data.frame(Ndiv)),as.data.frame(Ndiv))

  # 基于断面积ba
  sum.ba = tapply(data.Alive$ba,INDEX=list(data.Alive$Plot, data.Alive$SP), FUN=sum)
  sum.ba[is.na(sum.ba)]=0
  BAdiv = speciesdiv(sum.ba)
  BAdiv_out = data.frame("Plot"=rownames(as.data.frame(BAdiv)), as.data.frame(BAdiv)[,2:4])

  Species_result = left_join(Ndiv_out, BAdiv_out, by = 'Plot')
  colnames(Species_result) = c('Plot', 'Richness',
                               'N_GiniSimpson_Species_Div',	'N_Simpson_Species_Div',	'N_Shannon_Species_Div',	'Margalef_Div',	'Menhinick_Div',	'McIntosh_Div',
                               'BA_GiniSimpson_Species_Div',	'BA_Simpson_Species_Div',	'BA_Shannon_Species_Div')
  ####
  # 计算Chao1指数, 何潇添加，2024-4-7
  #  Chao1 = Richness+(F1(F1-1))/(2(F2+1))
  #  F1是只含有1个个体的树种数，F2是只含有2个个体的树种数，F1和F2相当于增加了珍贵树种的贡献
  F.temp = tapply(data.Alive$d,INDEX=list(data.Alive$Plot, data.Alive$SP), FUN=length)
  F.temp[is.na(F.temp)]=0
  Richness = apply(F.temp, 1, FUN = function(x) sum(x > 0))
  F1 = apply(F.temp, 1, FUN = function(x) sum(x == 1))
  F2 = apply(F.temp, 1, FUN = function(x) sum(x == 2))
  Chao1 = Richness + (F1*(F1-1))/(2*(F2+1))
  F.data = data.frame("Plot"=rownames(as.data.frame(F.temp)), Chao1)

  Species_result = dplyr::left_join(Species_result, F.data, by = 'Plot')

  #-------------------------------------------------------------------------------#
  # 径阶划分，此处并不智能
  # 小于5cm起测的样木均划分到一个径阶内
  # 大于51cm（对于2cm整化径阶）和53cm（对于4cm整化径阶）的样本均划分到一个径阶内
  # cut函数中 right = F,表示左闭右开
  # data.Alive$dbhclass2 <- cut(data.Alive$d, breaks = c(2*(4:27)-3,300), labels = c(2*(4:26)-2, '>51'), right = F)
  # data.Alive$dbhclass4 <- cut(data.Alive$d, breaks = c(4*(2:14)-3,300), labels = c(4*(2:13)-1, '>53'), right = F)

  # 径阶划分，周超凡-2022-8-15
  maxD = max(data.Alive$d, na.rm = T); minD = min(data.Alive$d, na.rm = T)
  maxDc2 = ceiling((maxD - minD + 0.001)/2)  #+0.001解决等于径级边界值的情况
  maxDc4 = ceiling((maxD - minD + 0.001)/4)  #+0.001解决等于径级边界值的情况
  data.Alive$dbhclass2 <- cut(data.Alive$d, breaks = c(2*(0:maxDc2)+minD), right = F)
  data.Alive$dbhclass4 <- cut(data.Alive$d, breaks = c(4*(0:maxDc4)+minD), right = F)
  #-------------------------------------------------------------------------------#

  ####
  ## 大小多样性 Size
  ## adiv包：speciesdiv()中输出的 GiniSimpson	Simpson	Shannon 3个指数计算结果
  ## 径阶划分包括2cm或4cm
  ## 可基于株数n或断面积ba计算
  ####
  # dbhclass2
  # 基于株数n
  counts.sp.dbhclass2 = tapply(data.Alive$d,INDEX=list(data.Alive$Plot, data.Alive$dbhclass2), FUN=length)
  counts.sp.dbhclass2[is.na(counts.sp.dbhclass2)]=0
  Ndiv_dbhclass2 = speciesdiv(counts.sp.dbhclass2)
  Ndiv_dbhclass2_out = data.frame("Plot"=rownames(as.data.frame(Ndiv_dbhclass2)), as.data.frame(Ndiv_dbhclass2)[,2:4])
  # 基于断面积ba
  sum.ba.dbhclass2 = tapply(data.Alive$ba,INDEX=list(data.Alive$Plot, data.Alive$dbhclass2), FUN=sum)
  sum.ba.dbhclass2[is.na(sum.ba.dbhclass2)]=0
  BAdiv_dbhclass2 = speciesdiv(sum.ba.dbhclass2)
  BAdiv_dbhclass2_out = data.frame("Plot"=rownames(as.data.frame(BAdiv_dbhclass2)),as.data.frame(BAdiv_dbhclass2)[,2:4])

  Size_dClass2_result = dplyr::left_join(Ndiv_dbhclass2_out, BAdiv_dbhclass2_out, by = 'Plot')
  colnames(Size_dClass2_result) = c('Plot',
                                    'N_GiniSimpson_Size2_Div',	'N_Simpson_Size2_Div',	'N_Shannon_Size2_Div',
                                    'BA_GiniSimpson_Size2_Div', 'BA_Simpson_Size2_Div',	'BA_Shannon_Size2_Div')
  # dbhclass4
  # 基于株数n
  counts.sp.dbhclass4 = tapply(data.Alive$d,INDEX=list(data.Alive$Plot, data.Alive$dbhclass4), FUN=length)
  counts.sp.dbhclass4[is.na(counts.sp.dbhclass4)]=0
  Ndiv_dbhclass4 = speciesdiv(counts.sp.dbhclass4)
  Ndiv_dbhclass4_out = data.frame("Plot"=rownames(as.data.frame(Ndiv_dbhclass4)), as.data.frame(Ndiv_dbhclass4)[,2:4])
  # 基于断面积ba
  sum.ba.dbhclass4 = tapply(data.Alive$ba,INDEX=list(data.Alive$Plot, data.Alive$dbhclass4), FUN=sum)
  sum.ba.dbhclass4[is.na(sum.ba.dbhclass4)]=0
  BAdiv_dbhclass4 = speciesdiv(sum.ba.dbhclass4)
  BAdiv_dbhclass4_out = data.frame("Plot"=rownames(as.data.frame(BAdiv_dbhclass4)),as.data.frame(BAdiv_dbhclass4)[,2:4])

  Size_dClass4_result = dplyr::left_join(Ndiv_dbhclass4_out, BAdiv_dbhclass4_out, by = 'Plot')
  colnames(Size_dClass4_result) = c('Plot',
                                    'N_GiniSimpson_Size4_Div',	'N_Simpson_Size4_Div',	'N_Shannon_Size4_Div',
                                    'BA_GiniSimpson_Size4_Div', 'BA_Simpson_Size4_Div',	'BA_Shannon_Size4_Div')

  ####
  ## 物种大小综合多样性 SS
  ## adiv包：speciesdiv()中输出的 GiniSimpson	Simpson	Shannon 3个指数计算结果
  ## 径阶划分包括2cm或4cm
  ## 可基于株数n或断面积ba计算
  ####
  data.Alive$SP_dbhclass2 = paste0(data.Alive$SP, data.Alive$dbhclass2)
  data.Alive$SP_dbhclass4 = paste0(data.Alive$SP, data.Alive$dbhclass4)

  # SP_dbhclass2
  # 基于株数n
  counts.sp.sp.dbhclass2 = tapply(data.Alive$d,INDEX=list(data.Alive$Plot, data.Alive$SP_dbhclass2), FUN=length)
  counts.sp.sp.dbhclass2[is.na(counts.sp.sp.dbhclass2)]=0
  Ndiv_sp_dbhclass2 = speciesdiv(counts.sp.sp.dbhclass2)
  Ndiv_sp_dbhclass2_out = data.frame("Plot"=rownames(as.data.frame(Ndiv_sp_dbhclass2)), as.data.frame(Ndiv_sp_dbhclass2)[,2:4])
  # 基于断面积ba
  sum.ba.sp.dbhclass2 = tapply(data.Alive$ba,INDEX=list(data.Alive$Plot, data.Alive$SP_dbhclass2), FUN=sum)
  sum.ba.sp.dbhclass2[is.na(sum.ba.sp.dbhclass2)]=0
  BAdiv_SP_dbhclass2 = speciesdiv(sum.ba.sp.dbhclass2)
  BAdiv_SP_dbhclass2_out = data.frame("Plot"=rownames(as.data.frame(BAdiv_SP_dbhclass2)),as.data.frame(BAdiv_SP_dbhclass2)[,2:4])

  SS_dClass2_result = dplyr::left_join(Ndiv_sp_dbhclass2_out, BAdiv_SP_dbhclass2_out, by = 'Plot')
  colnames(SS_dClass2_result) = c('Plot',
                                  'N_GiniSimpson_SS2_Div',	'N_Simpson_SS2_Div',	'N_Shannon_SS2_Div',
                                  'BA_GiniSimpson_SS2_Div', 'BA_Simpson_SS2_Div',	'BA_Shannon_SS2_Div')
  # SP_dbhclass4
  # 基于株数n
  counts.sp.sp.dbhclass4 = tapply(data.Alive$d,INDEX=list(data.Alive$Plot, data.Alive$SP_dbhclass4), FUN=length)
  counts.sp.sp.dbhclass4[is.na(counts.sp.sp.dbhclass4)]=0
  Ndiv_sp_dbhclass4 = speciesdiv(counts.sp.sp.dbhclass4)
  Ndiv_sp_dbhclass4_out = data.frame("Plot"=rownames(as.data.frame(Ndiv_sp_dbhclass4)), as.data.frame(Ndiv_sp_dbhclass4)[,2:4])
  # 基于断面积ba
  sum.ba.sp.dbhclass4 = tapply(data.Alive$ba,INDEX=list(data.Alive$Plot, data.Alive$SP_dbhclass4), FUN=sum)
  sum.ba.sp.dbhclass4[is.na(sum.ba.sp.dbhclass4)]=0
  BAdiv_SP_dbhclass4 = speciesdiv(sum.ba.sp.dbhclass4)
  BAdiv_SP_dbhclass4_out = data.frame("Plot"=rownames(as.data.frame(BAdiv_SP_dbhclass4)),as.data.frame(BAdiv_SP_dbhclass4)[,2:4])

  SS_dClass4_result = dplyr::left_join(Ndiv_sp_dbhclass4_out, BAdiv_SP_dbhclass4_out, by = 'Plot')
  colnames(SS_dClass4_result) = c('Plot',
                                  'N_GiniSimpson_SS4_Div',	'N_Simpson_SS4_Div',	'N_Shannon_SS4_Div',
                                  'BA_GiniSimpson_SS4_Div', 'BA_Simpson_SS4_Div',	'BA_Shannon_SS4_Div')

  ## 控制输出
   if (any(!Index%in% c('Species', 'Size', 'SS', 'ALL_2', 'ALL_4', 'ALL')))
      stop("Your choice for Index is not available")
  else if(Index=='Species'){
    return(Species_result)
  }
  else if(Index=='Size'){
    if(!is.null(dClass)){
      if(dClass==2){
        Size_result = Size_dClass2_result
      }
      else if(dClass==4){
        Size_result = Size_dClass4_result
      }
    }
    else{
      Size_result = dplyr::left_join(Size_dClass2_result, Size_dClass4_result, by = 'Plot')
    }
    return(Size_result)
  }
  else if(Index=='SS'){
    if(!is.null(dClass)){
      if(dClass==2){
        SS_result = SS_dClass2_result
      }
      else if(dClass==4){
        SS_result = SS_dClass4_result
      }
    }
    else{
      SS_result = dplyr::left_join(SS_dClass2_result, SS_dClass4_result, by = 'Plot')
    }
    return(SS_result)
  }
  else if(Index=='ALL_2'){
    result2 = dplyr::left_join(Species_result, Size_dClass2_result, by = 'Plot')
    result2 = dplyr::left_join(result2, SS_dClass2_result, by = 'Plot')
    return(result2)
  }
  else if(Index=='ALL_4'){
    result4 = dplyr::left_join(Species_result, Size_dClass4_result, by = 'Plot')
    result4 = dplyr::left_join(result4, SS_dClass4_result, by = 'Plot')
    return(result4)
  }
  else if(Index=='ALL'){
    Size_result = dplyr::left_join(Size_dClass2_result, Size_dClass4_result, by = 'Plot')
    SS_result = dplyr::left_join(SS_dClass2_result, SS_dClass4_result, by = 'Plot')
    result = dplyr::left_join(Species_result, Size_result, by = 'Plot')
    result = dplyr::left_join(result, SS_result, by = 'Plot')
    return(result)
  }
}
