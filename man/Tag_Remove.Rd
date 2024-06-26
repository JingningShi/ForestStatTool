% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Tag_Remove.R
\name{Tag_Remove}
\alias{Tag_Remove}
\title{重复标签处理2}
\usage{
Tag_Remove(Data, Plot, Tag, D)
}
\arguments{
\item{Data}{（必须）数据框。存在标签重复的样地林木数据框。}

\item{Plot}{（必须）向量，建议为字符型。样地号。}

\item{Tag}{（必须）向量，建议为字符型。样地林木标签。}

\item{D}{（必须）数字向量。样地林木胸径。单位：cm。}
}
\value{
1）输出结果为一个数据框。

包括输入的Data数据框的所有列，以及函数可能产生的新列。

2）新输出变量名：

newTag：重新编号的新标签，格式为“原标签-编号”
}
\description{
Tag_Remove通过删除标签重复的较小胸径的林木而保留标签重复的最大胸径的林木来解决标签重复问题。
}
\details{
1）当各样地均不存在重复标签时，会输出“There are no duplicated tags!”的提醒，此时将不会再输出新数据框。

2）可对多个样地组成的数据框进行重复标签批量处理。
}
\examples{
## 加载内置数据
data(ForestStatTool)
##rawdata是一个存在重复标签的样地原始数据框
rawdata
##Tag_Remove()通过保留重复标签中胸径最大的林木来解决标签重复问题。
##参数D有空值，删除存在空值的行
rawdata0 <- subset(rawdata,!is.na(rawdata$D))
newdata2 <- Tag_Remove(Data = rawdata0, Plot = rawdata0$Plot, Tag = rawdata0$Tag, D = rawdata0$D)
newdata2
##检查newdata2中的Tag是否还存在重复标签
Tag_Remove(Data = newdata2, Plot = newdata2$Plot, Tag = newdata2$Tag, D = newdata2$D)
#输出“There are no duplicated tags!”提醒，说明所有样地中的样木都不存在重复标签
}
\author{
Chaofan Zhou: cfzhou2021@163.com

Xiao He: hexiao@ifrit.ac.cn

Guangshuang Duan: oliverdgs@163.com

Jingning Shi: shijingning@hebau.edu.cn

Xiangdong Lei: xdlei@ifrit.ac.cn
}
