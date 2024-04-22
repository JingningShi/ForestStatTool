#' @title Distance-independent competition indices. 与距离无关的竞争指数_9个

#' @author Guangshuang Duan: oliverdgs@163.com
#' @author Xiao He: hexiao@ifrit.ac.cn
#' @author Chaofan Zhou: cfzhou2021@163.com
#' @author Mengli Zhou: 1257304343@qq.com
#' @author Jingning Shi: shijingning@hebau.edu.cn
#' @author Xiangdong Lei: xdlei@ifrit.ac.cn

#' @usage CI.nondist(Data = NULL, Plot, Tag, D, SP = NULL, H = NULL, HCB = NULL, CR = NULL, S = 0.06, Bind = FALSE)

#' @description 此函数实现了利用单木测量因子计算与距离无关的单木竞争指数，包括与距离无关的9个变量。

#' @details 1）输入数据可以选择两种方法：第一种是按照参数Plot, Tag, D, SP, H, HCB和CR逐个输入，其中SP, H, HCB和CR可选择不输入；第二种是首先修改数据的列名，包括Plot, Tag, D, SP, H, CR，其中SP, H, HCB和CR可不包含在数据中，然后通过Data参数一次性输入。
#' @details 2）如果SP不输入则只计算总竞争，输入SP后会同时计算9个竞争指数的总竞争和种内种间竞争。
#' @details 3）如果H不输入则CI8不计算，如果H、HCB和CW中有一个不输入则CI9不计算。
#' @details 4）与距离无关的竞争指数如下：
#' @details C1	大于对象木胸径的其他林木的胸高断面积之和
#' @details C2	胸径不等于对象木胸径的其他林木的胸径之和与对象木胸径的比
#' @details C3	对象木胸径与样地中林木的平均胸径（均方直径）的比
#' @details C4	对象木胸径与样地中林木的最大胸径的比
#' @details C5	对象木的胸高断面积与样地中林木的平均胸高断面积的比
#' @details C6	对象木的胸高断面积与样地中林木的最大胸高断面积的比
#' @details C7	样地中所有胸径不等于对象木胸径的其他林木的胸高断面积之和与对象木胸高断面积的比
#' @details C8	样地中大于对象木胸径的其他林木的胸高断面积之和与样地单位面积胸高断面积和相对空间指数的比。注：需要单木树高
#' @details C9	单位面积的样地中活立木最大树冠高度为对象木树高的66%的竞争木的冠幅面积之和。注：需要单木树高、枝下高、冠幅

#' @param Data （可选）数据框，用于存放单木数据。其变量必须含有Plot、Tag、D。在Bind = TRUE 时，Data必须输入。
#' @param Plot （可选）向量，样地号。可以是数值型，也可以是字符型。
#' @param Tag （可选）向量，林木标签（编号）。可以是数值型，也可以是字符型。
#' @param D （必须）数字向量。单木胸径，单位：cm。
#' @param SP （可选）向量，可以是数字型，也可以是字符型。树种。
#' @param H （可选）向量。样地林木胸径。单位：m。默认为空。
#' @param HCB （可选）向量。单木枝下高。单位：m。默认为空。
#' @param CR （可选）向量。单木冠幅半径。单位：m。默认为空。
#' @param S （必须）数字。样地面积。单位：ha。默认为0.06ha。
#' @param Bind （必须）逻辑向量。是否将输出的空间结构指标与Data进行匹配合并。默认为FALSE，不进行合并。当Bind = TRUE时，将空间结构指标按照样地号和林木标签与Data进行匹配后合并输出。

#' @return 1）输出结果为一个数据框。如果Bind=FLASE，数据框仅包含样地号Plot、林木标签Tag和竞争指数；若Bind=TRUE且Data不为空，数据框包含Data的所有林木信息和竞争指数。
#' @return 2）输出变量名：
#' @return Plot：样地号
#' @return Tag：林木标签
#' @return CI1~CI9：9个与距离无关的总竞争指数
#' @return CI1_intra~ CI9_intra：9个与距离无关的种内竞争指数
#' @return CI1_inter~CI9_inter：9个与距离无关的种间竞争指数

#' @name CI.nondist
#' @import stats

#' @references Wykoff, W.R., Crookston, N.L., Stage, A.R. (1982). User’s guide to the stand prognosis model (Vol. 133). United States Department of Agriculture, Forest Service, Intermountain Forest and Range Experiment Station. General Technical Report. INT-133, pp. 119.
#' @references Lorimer, C.G. (1983). Tests of age-independent competition indices for individual trees in natural hardwood stands. Forest Ecology and Management, 6: 343-360.
#' @references Hamilton, D. A. (1986). A logistic model of mortality in thinned and unthinned mixed conifer stands of northern Idaho. Forest Science, 32(4), 989-1000.
#' @references Tomé, M., Burkhart, H. E. (1989). Distance-dependent competition measures for predicting growth of individual trees. Forest Science, 35(3), 816-831.
#' @references Corona, P., Ferrara, A. (1989). Individual competition indices for conifer plantations. Agriculture Ecosystems and Environment, 27(1-4), 429-437.
#' @references Schröder, J., Gadow, K.V. (1999). Testing a new competition index for Maritime pine in northwestern Spain. Canandian Journal of Forest Research, 29, 280-283.
#' @references Nagel, J. (1999). Konzeptionelle Überlegungen zum schrittweisen Aufbau eineswaldwachstumskundlichen Simulationssystems für Nordwestdeutschland. Schriften aus der Forstlichen Fakultät der Universität Göttingen und der Niedersächsischen Forstlichen Versuchsanstalt, Band 128, J.D. Sauerländer's Verlag, Frankfurt am Main, pp. 122.

#' @export CI.nondist

#' @importFrom magrittr "%>%"
#' @importFrom dplyr arrange filter left_join mutate select
#' @importFrom plyr create_progress_bar progress_text

#' @examples ## 加载内置数据
#' @examples data(ForestStatTool)
#' @examples #bigplot是一个100m×100m的矩形样地的林木信息，使用参数逐个输入林木信息，计算与距离无关的竞争指数
#' @examples a1 <- CI.nondist (Plot = bigplot$Plot, Tag = bigplot$Tag, D = bigplot$D, SP = bigplot$SP, S = 1)
#' @examples head(a1)
#' @examples #使用Data参数直接输入，由于Data中除了必须输入的变量外，还含有SP、H、HCB和CR,会计算更多竞争指数
#' @examples a2 <- CI.nondist (Data = bigplot, S=1)
#' @examples head(a2)
#' @examples #cdf是一个半径13.82m的圆形样地的林木信息，使用参数逐个输入林木信息，计算与距离无关的竞争指数
#' @examples a3 <- CI.nondist (Plot = cdf$Plot, Tag = cdf$Tag, D = cdf$D, SP = cdf$SP, S = 0.06)
#' @examples head(a3)


CI.nondist <- function(Data=NULL, Plot, Tag, D, SP = NULL, H = NULL, HCB = NULL, CR = NULL, S=0.06, Bind = FALSE){
 # if(!require("plyr")){
 #   install.packages("plyr")
 #   library(plyr)
 # }
 # if(!require("dplyr")){
 #   install.packages("dplyr")
 #   library(dplyr)
 # }

  #---------计算优势高、优势径,周梦丽论文推荐的方法
  UHUD.f <- function(data,S){
    data<-arrange(data,desc(D))
    n1<-floor(1.6*S*100-0.6)
    n2<-floor(1.6*S*100-0.6)+1
    H1<-mean(data$H[1:n1], na.rm = T)
    D1<-mean(data$D[1:n1], na.rm = T)
    H2<-mean(data$H[1:n2], na.rm = T)
    D2<-mean(data$D[1:n2], na.rm = T)
    UD<-D1+(D2-D1)/(n2-n1)*0.4
    UH<-H1+(H2-H1)/(n2-n1)*0.4
    UHUD<-c(UD,UH)
    return(UHUD)
  }
  #---------求和
  sum.f <- function(data){
    if(class(data)=="data.frame"){
      if(nrow(data)==0) {
        sum<-0
      } else {
        sum<-sum(data, na.rm = T)
      }
    } else {
      if(length(data)==0) {
        sum<-0
      } else {
        sum<-sum(data, na.rm = T)
      }
    }
    return(sum)
  }
  #---------求面积
  area.f <- function(data){
    if(class(data)=="data.frame"){
      if(nrow(data)==0) {
        area<-0
      } else {
        area<-pi*sum(data^2, na.rm = T)
      }
    } else {
      if(length(data)==0) {
        area<-0
      } else {
        area<-pi*sum(data^2, na.rm = T)
      }
    }
    return(area)
  }

  #---------作用对象是一个样地,data需要字段胸径D、树高H
  #---------分有无树高两种情况
  AddIndex.f <- function(data,S){
    if ("H" %in% colnames(data)){
      temp <- UHUD.f(data,S=S)[2]
      data <- mutate(data,S=S, N=nrow(data), Hdom=temp,
                            Dg=sqrt(mean(data$D^2, na.rm = T)), Dmax=max(data$D, na.rm = T),
                            BA=pi/40000*data$D^2, BAmean=mean(BA, na.rm = T), BAmax=max(BA, na.rm = T), BAsum=sum(BA, na.rm = T)/S,
                            RS=sqrt(10000*S/N)/Hdom)
    } else {
      data <- mutate(data,S=S, N=nrow(data),
                            Dg=sqrt(mean(data$D^2, na.rm = T)), Dmax=max(data$D, na.rm = T),
                            BA=pi/40000*data$D^2, BAmean=mean(BA, na.rm = T), BAmax=max(BA, na.rm = T), BAsum=sum(BA, na.rm = T)/S)
    }
    return(data)
  }

  # 统计NA的比例，何潇-2022-12-8
  Na.f=function(x){ round(sum(is.na(x))/length(x)*100, 2) }
  #数据输入可选择数据框，也可以是向量
  # 补充没有H不计算CI8，何潇-2022-11-27
  # 数据组织统一，何潇2022-12-10
  if(is.null(Data)){
    data <- data.frame("Plot"=Plot,"Tag"=Tag,"D"=D)
    if(!is.null(SP)){
      data$SP <- SP
    }
    if(!is.null(H)){
      data$H <- H
    }
    if(!is.null(HCB)){
      data$HCB <- HCB
    }
    if(!is.null(CR)){
      data$CR <- CR
    }
    if(Bind==TRUE){
      stop("Missing 'Data', 'Bind = TRUE' is not applicable.")
    }
  } else{ # 以Data的形式输入时，要进行变量名检查，何潇-200-12-8
    if(!"Plot"%in%colnames(Data)){
      stop("Missing column name of 'Plot' in Data.")
    }
    if(!"Tag"%in%colnames(Data)){
      stop("Missing column name of 'Tag' in Data.")
    }
    if(!"D"%in%colnames(Data)){
      stop("Missing column name of 'D' in Data.")
    }
    data <- select(Data, Plot, Tag, D) # Data中必须含有的变量
    # 竞争指数中是否计算种内种间根据是否输入SP来判断，何潇2022-12-10
    if("SP"%in%colnames(Data)){
      data$SP <- Data$SP   ## 发现bug,数据组织有问题，何潇修改2024-3-9
    }
    # 判断是否输入H、HCB和CR变量，何潇2022-12-10
    if("H"%in%colnames(Data)){
      data$H <- Data$H
    }
    if("HCB"%in%colnames(Data)){
      data$HCB <- Data$HCB
    }
    if("CR"%in%colnames(Data)){
      data$CR <- Data$CR
    }
  }

  # 以下3个条件判断有误，何潇修改输出条件，2024-3-9
  if((!"SP"%in%colnames(data))|nlevels(factor(data$SP))==1){
    cat("The 'CIx_intra' and 'CIx_inter' not calculate because 'SP' don't input or the levels of 'SP' equal to 1.\n")
  }
  if(!"H"%in%colnames(data)){
    cat("The 'CI8' not calculate because 'H' don't input.\n")
  }
  if(!all("H"%in%colnames(data), "HCB"%in%colnames(data), "CR"%in%colnames(data))){
    cat("The 'CI9' not calculate because 'H', 'HCB' or 'CR' don't input.\n")
  }

  #主要变量合理性检测
  # 修改为变量有缺失值为打印，而不是停止，下面有删除NA的操作。何潇2022-12-8
  # 如果参数的类型不对，如D不是数值型，则停止计算。何潇2022-12-8
  if(sum(is.na(data$Plot))!=0){
    cat(paste0(Na.f(data$Plot), "% Missing value (or NA) in 'Plot'.\n"))
  }
  if(sum(is.na(data$Tag))!=0){
    cat(paste0(Na.f(data$Tag), "% Missing value (or NA) in 'Tag'.\n"))
  }
  ##
  if(!is.numeric(data$D)){
    stop("'D' must be numeric.")
  }else if(sum(is.na(data$D))!=0){
    cat(paste0(Na.f(data$D), "% Missing value (or NA) in 'D'.\n"))
  }
  ##
  # SP只检查是否有空值
  if(!is.null(data$SP)){
    if(T %in% is.na(data$SP)){
      cat(paste0(Na.f(data$SP), "Missing value (or NA) in 'SP'.\n"))
    }
  }
  # H、HCB、CR只有输入了才检查
  if(!is.null(data$H)){
    if(!is.numeric(data$H)){
      stop("'H' must be numeric.")
    }else if(T %in% is.na(data$H)){
      cat(paste0(Na.f(data$H), "% Missing value (or NA) in 'H'.\n"))
    }
  }
  if(!is.null(data$HCB)){
    if(!is.numeric(data$HCB)){
      stop("'HCB' must be numeric.")
    }else if(T %in% is.na(data$HCB)){
      cat(paste0(Na.f(data$HCB), "% Missing value (or NA) in 'HCB'.\n"))
    }
  }
  if(!is.null(data$CR)){
    if(!is.numeric(data$CR)){
      stop("'CR' must be numeric.")
    }else if(T %in% is.na(data$CR)){
      cat(paste0(Na.f(data$HCB), "% Missing value (or NA) in 'CR'.\n"))
    }
  }

  #删除NA值防止无法运算，何潇-2022-11-27
  data <- na.omit(data)
  if(nrow(data)==0){
    temp = apply(data, 2, function(x) sum(!is.na(x))==0)
    tempName = names(which(temp == TRUE))
    if(length(tempName)==1){
      stop(paste0("The column of", tempName, "is NA."))
    }else{
      stop(paste0("The columns of", tempName, "are NA."))
    }
  }

  # 总树种个数
  N_sp <- nlevels(factor(data$SP))
  #以Plot分组计算竞争指标
  N <- nlevels(factor(data$Plot))
  # 增加进度条，何潇-2022-11-28
  cat("Start calculating: \n")
  progress.bar <- create_progress_bar("text")  #plyr包中的create_progress_bar函数创建一个进度条
  progress.bar$init(N)   #设置任务数，几个样地

  # 计算总竞争
  data_jieguo<-data.frame()
  data_jieguo2<-data.frame() #用于记录种内和种间竞争的结果
  for (j in 1:N) {#先循环样地后判断，周超凡2024-3-4
    data1 <- filter(data, Plot==levels(factor(Plot))[j]) # 所有都要计算总竞争，提取出来，减少多余代码，周超凡2024-3-4
    #重复标签检查, 周超凡2022-12-11
    # if(TRUE %in% duplicated(data1$Tag)) {} # 何潇修改判断语句写法,2024-4-10
    if( any(duplicated(data1$Tag)) ){
      stop("There are tags duplicated in data, please use the function 'Tag_Divide' or 'Tag_Remove' to solve, and use the new column of 'newTag' to analyse.")
    }
    #向量化运算
    data1 <- AddIndex.f(data1, S=S)
    data1 <- mutate(data1, CI3=D/Dg, CI4=D/Dmax, CI5=BA/BAmean, CI6=BA/BAmax)
    #无法向量化的采用循环
    for (i in 1:nrow(data1)){
      data1$CI1[i]<-select(filter(data1,BA > BA[i]),BA) %>% sum.f()/S
      data1$CI2[i]<-select(filter(data1,D != D[i]),D) %>% sum.f()/data1$D[i]
      data1$CI7[i]<-select(filter(data1,D != D[i]),BA) %>% sum.f()/data1$BA[i]
      #判断是否计算CI8
      if(!is.null(data1$H)){
        data1$CI8[i]<-select(filter(data1, D > D[i]), BA) %>% sum.f()/data1$BAsum[i]/data1$RS[i]
        #判断是否计算CI9
        if(!is.null(data1$HCB)&!is.null(data1$CR)){
          data1$CI9[i]<-select(filter(data1, H-HCB <= 0.66*H[i]), CR) %>% area.f()/S
        } else {
          data1$CI9[i]<-NA
        }
      } else {
        data1$CI8[i]<-NA
        data1$CI9[i]<-NA
      }
    }
    data_jieguo<-rbind(data_jieguo,data1)
    #计算种内竞争、种间竞争
    if(!(is.null(data$SP)|(N_sp==1))){
      # 计算种内种间竞争
      N_sp_j <- nlevels(factor(data1$SP))
      for (k in 1:N_sp_j) {
        # 向量化计算
        # data2 保存相同树种的数据
        data2 <- filter(data1, SP==levels(factor(SP))[k])
        data2 <- AddIndex.f(data2,S=S)
        # data3 保存不同树种的数据
        data3 <- filter(data1, SP!=levels(factor(SP))[k])
        data3 <- AddIndex.f(data3,S=S)
        # 向量化计算结果
        data2 <- mutate(data2,
                        CI3_intra=D/Dg, CI4_intra=D/Dmax, CI5_intra=BA/BAmean, CI6_intra=BA/BAmax,
                        CI3_inter=D/mean(data3$Dg), CI4_inter=D/mean(data3$Dmax), CI5_inter=BA/mean(data3$BAmean), CI6_inter=BA/mean(data3$BAmax))
        # 非向量化计算
        for (i in 1:nrow(data2)){
          data2$CI1_intra[i]<-select(filter(data2, BA > BA[i]), BA) %>% sum.f()/S
          data2$CI1_inter[i]<-select(filter(data3, BA > data2$BA[i]), BA) %>% sum.f()/S
          data2$CI2_intra[i]<-select(filter(data2, D != D[i]), D) %>% sum.f()/data2$D[i]
          data2$CI2_inter[i]<-select(filter(data3, D != data2$D[i]), D) %>% sum.f()/data2$D[i]
          data2$CI7_intra[i]<-select(filter(data2, D != D[i]), BA) %>% sum.f()/data2$BA[i]
          data2$CI7_inter[i]<-select(filter(data3, D != data2$D[i]), BA) %>% sum.f()/data2$BA[i]
          #判断是否计算CI8_intra
          if(!is.null(data2$H)){
            # 将data2和data3中的优势高Hdom和相对空间指数RS和林分断面积转换为data1中样地尺度计算的结果
            data2 <- data2 %>% select(-c(N, Hdom, RS, BAsum)) %>%
              left_join(data1 %>% select(Plot, Tag, N, Hdom, RS, BAsum), by=c("Plot", "Tag"))
            data3 <- data3 %>% select(-c(N, Hdom, RS, BAsum)) %>%
              left_join(data1 %>% select(Plot, Tag, N, Hdom, RS, BAsum), by=c("Plot", "Tag"))
            # 开始计算
            data2$CI8_intra[i]<-select(filter(data2, D > D[i]),BA) %>% sum.f()/data2$BAsum[i]/data2$RS[i]
            data2$CI8_inter[i]<-select(filter(data3, D > data2$D[i]),BA) %>% sum.f()/data2$BAsum[i]/data2$RS[i]
            #判断是否计算种内CI9_intra
            if(!is.null(data2$HCB)&!is.null(data2$CR)){
              data2$CI9_intra[i]<-select(filter(data2, H-HCB <= 0.66*H[i]), CR) %>% area.f()/S
              data2$CI9_inter[i]<-select(filter(data3, H-HCB <= 0.66*data2$H[i]), CR) %>% area.f()/S
            } else {
              data2$CI9_intra[i]<-NA
              data2$CI9_inter[i]<-NA
            }
          } else {
            data2$CI8_intra[i]<-NA
            data2$CI8_inter[i]<-NA
            data2$CI9_intra[i]<-NA
            data2$CI9_inter[i]<-NA
          }
        }
        # 结果合并，段光爽修改，2022.12.01
        data_jieguo2<-rbind(data_jieguo2, data2)
      }
    }
    # 每计算完一个样地打印一次结果
    cat("\n")
    cat(paste("Calculated" ,j, ": Plot =", levels(factor(data_jieguo$Plot))[j]), "\n")
    progress.bar$step() #输出进度条
  }
  cat("\n")
  data_jieguo <- select(data_jieguo,Plot,Tag, paste("CI",1:9,sep = ""))


  # 根据输入的变量判断总竞争data_jieguo的输出结果，何潇2022-12-2
  # 判断条件错误，应该是判断data中是否存在变量名，何潇2024-3-9
  if(all(c("H") %in% names(data)) & !all(c("HCB","CR") %in% names(data))){
    data_jieguo <- select(data_jieguo,Plot,Tag, paste("CI",1:8,sep = ""))
  }else {#else if后有else导致出错，直接调用else，修改成重新判断，周超凡2024-3-4
    if(all(c("H", "HCB","CR") %in% names(data))){
      data_jieguo <- select(data_jieguo,Plot,Tag, paste("CI",1:9,sep = ""))
    }else{
      data_jieguo <- select(data_jieguo,Plot,Tag, paste("CI",1:7,sep = ""))
    }
  }


  # 当不输入树种SP或输入的SP中只有一个树种时，data_jieguo2不存在，何潇2022-12-01
  # 根据输入的变量判断种内种间竞争data_jieguo2输出结果，何潇2022-12-2
  # 判断条件错误，应该是判断data中是否存在变量名，何潇2024-3-9
  if(!is.null(data$SP)|(N_sp==1)){
    if(all(c("H") %in% names(data)) & !all(c("HCB","CR") %in% names(data))){
      data_jieguo2<-data_jieguo2 %>%
        select(Plot,Tag, paste("CI",1:8,"_intra", sep = ""), paste("CI",1:8,"_inter", sep = ""))
    }else{#else if后有有else导致出错，直接调用else，修改成重新判断，周超凡2024-3-4
      if(all(c("H", "HCB","CR") %in% names(data))){
        data_jieguo2<-data_jieguo2 %>%
          select(Plot,Tag, paste("CI",1:9,"_intra", sep = ""), paste("CI",1:9,"_inter", sep = ""))
      }else{
        data_jieguo2<-data_jieguo2 %>%
          select(Plot,Tag, paste("CI",1:7,"_intra", sep = ""), paste("CI",1:7,"_inter", sep = ""))
      }
    }
    data_jieguo <- data_jieguo %>% left_join(data_jieguo2, by=c("Plot","Tag"))
  }
  #判断是否与原数据框合并
  if(Bind&!is.null(Data)){
    data_jieguo <- left_join(Data, data_jieguo, by = c("Plot","Tag"))
  }
  return(data_jieguo)
}

