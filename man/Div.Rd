% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Div.R
\name{Div}
\alias{Div}
\title{多样性计算}
\usage{
Div(Plot, D, SP, TreeType = NULL, Dmin = 5,  dClass = NULL, Index = 'Species')
}
\arguments{
\item{Plot}{（必须）向量，可以是数值型，也可以是字符型。样地号。}

\item{D}{（必须）数值向量。样地林木胸径。单位：cm。}

\item{SP}{（可选）向量，可以是数值型，也可以是字符型。树种。}

\item{TreeType}{（可选）数值向量。检尺类型。12和14表示枯死木，13和50表示采伐木，其余值均认为是活立木，可输入全部检尺类型，程序将自动计算活立木（检尺类型不等于12、13、14和50）的生物多样性因子。默认为空，即输入的数据均为活立木。}

\item{Dmin}{数值。起测胸径，默认为5cm。}

\item{dClass}{数值。径阶，2或4。其值为空时会同时输出2cm和4cm径阶划分所计算的结构多样性指数，默认为空。}

\item{Index}{字符串。指标名称，取值为'Species'、'Size'、'SS'、'ALL_2'、'ALL_4'、'ALL'之一。'Species'，计算物种多样性/均匀度指数；'Size'，计算大小多样性/均匀度指数；'SS'，计算物种大小综合多样性/均匀度指数；搭配dClass参数可实现只计算2cm或只计算4cm径阶划分的大小多样性/均匀度指数、综合多样性/均匀度指数；'ALL_2'，计算物种多样性/均匀度指数、2cm径阶划分的大小多样性/均匀度指数以及2cm径阶划分的综合多样性/均匀度指数；'ALL_4'，计算物种多样性/均匀度指数、4cm径阶划分的大小多样性/均匀度指数以及4cm径阶划分的综合多样性/均匀度指数；'ALL'，计算所有多样性/均匀度指数。当Index取值为'ALL_2'、'ALL_4'、 'ALL'之一时，dClass参数不起作用。默认为'Species'。}
}
\value{
输出结果为一个数据框。

2）输出变量名：

Plot：样地号

Richness：树种丰富度指数

N_GiniSimpson_Species_Div：Pi为株数占比时的GiniSimpson树种多样性指数

N_Simpson_Species_Div：Pi为株数占比时的Simpson树种多样性指数

N_Shannon_Species_Div：Pi为株数占比时的Shannon树种多样性指数

Margalef_Div：Margalef丰富度指数

Menhinick_Div：Menhinick丰富度指数

McIntosh_Div：McIntosh丰富度指数

BA_GiniSimpson_Species_Div：Pi为断面积占比时的GiniSimpson树种多样性指数

BA_Simpson_Species_Div：Pi为断面积占比时的Simpson树种多样性指数

BA_Shannon_Species_Div：Pi为断面积占比时的Shannon树种多样性指数

Chao1：Chao1指数

N_GiniSimpson_Size2_Div：Pi为株数占比时的GiniSimpson大小多样性指数（按2cm划分径阶）

N_Simpson_Size2_Div：Pi为株数占比时的Simpson大小多样性指数（按2cm划分径阶）

N_Shannon_Size2_Div：Pi为株数占比时的Shannon大小多样性指数（按2cm划分径阶）

BA_GiniSimpson_Size2_Div：Pi为断面积占比时的GiniSimpson大小多样性指数（按2cm划分径阶）

BA_Simpson_Size2_Div：Pi为断面积占比时的Simpson大小多样性指数（按2cm划分径阶）

BA_Shannon_Size2_Div：Pi为断面积占比时的Shannon大小多样性指数（按2cm划分径阶）

N_GiniSimpson_Size4_Div：Pi为株数占比时的GiniSimpson大小多样性指数（按4cm划分径阶）

N_Simpson_Size4_Div：Pi为株数占比时的Simpson大小多样性指数（按4cm划分径阶）

N_Shannon_Size4_Div：Pi为株数占比时的Shannon大小多样性指数（按4cm划分径阶）

BA_GiniSimpson_Size4_Div：Pi为断面积占比时的GiniSimpson大小多样性指数（按4cm划分径阶）

BA_Simpson_Size4_Div：Pi为断面积占比时的Simpson大小多样性指数（按4cm划分径阶）

BA_Shannon_Size4_Div：Pi为断面积占比时的Shannon大小多样性指数（按4cm划分径阶）

N_GiniSimpson_SS2_Div：Pi为株数占比时的GiniSimpson综合多样性指数（按2cm划分径阶）

N_Simpson_SS2_Div：Pi为株数占比时的Simpson综合多样性指数（按2cm划分径阶）

N_Shannon_SS2_Div：Pi为株数占比时的Shannon综合多样性指数（按2cm划分径阶）

BA_GiniSimpson_SS2_Div：Pi为断面积占比时的GiniSimpson综合多样性指数（按2cm划分径阶）

BA_Simpson_SS2_Div：Pi为断面积占比时的Simpson综合多样性指数（按2cm划分径阶）

BA_Shannon_SS2_Div：Pi为断面积占比时的Shannon综合多样性指数（按2cm划分径阶）

N_GiniSimpson_SS4_Div：Pi为株数占比时的GiniSimpson综合多样性指数（按4cm划分径阶）

N_Simpson_SS4_Div：Pi为株数占比时的Simpson综合多样性指数（按4cm划分径阶）

N_Shannon_SS4_Div：Pi为株数占比时的Shannon综合多样性指数（按4cm划分径阶）

BA_GiniSimpson_SS4_Div：Pi为断面积占比时的GiniSimpson综合多样性指数（按4cm划分径阶）

BA_Simpson_SS4_Div：Pi为断面积占比时的Simpson综合多样性指数（按4cm划分径阶）

BA_Shannon_SS4_Div：Pi为断面积占比时的Shannon综合多样性指数（按4cm划分径阶）
}
\description{
此函数实现了输入单木胸径和树种等调查变量，统计多样性指数。

具体还包括物种多样性、（非空间）结构多样性等。
}
\examples{
## 加载内置数据
data(ForestStatTool)
#bigplot是一个100m×100m的矩形样地的林木信息，首先利用Plot_Divide函数将bigplot划分为后的25个20m×20m的小样地的合并数据框。
#除原始变量外，还包含的新列有：分割后的新样地名subplot，小样地的新X坐标subX和小样地的新Y坐标subY
subplot1 <- Plot_Divide(Data = bigplot, Plot = bigplot$Plot, X = bigplot$X, Y = bigplot$Y,
                      Num_xy = c(5,5), Range_xy = c(100,100))
#计算不同小样地的多样性指标
a1 = Div(Plot = subplot1$subplot, D = subplot1$D, SP = subplot1$SP, dClass = NULL, Index = 'ALL')
head(a1)
}
\author{
Xiao He: hexiaonuist@163.com

Chaofan Zhou: cfzhou2021@163.com

Guangshuang Duan: oliverdgs@163.com

Jingning Shi: shijingning@hebau.edu.cn

Xiangdong Lei: xdlei@ifrit.ac.cn
}
