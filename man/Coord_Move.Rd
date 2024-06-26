% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Coord_Move.R
\name{Coord_Move}
\alias{Coord_Move}
\title{重复坐标点处理1}
\usage{
Coord_Move(Data, Plot, X, Y, Origin = c(0,0), Range_xy)
}
\arguments{
\item{Data}{（必须）数据框。存在坐标重复的样地林木数据框。}

\item{Plot}{（必须）向量，建议为字符型。样地号。}

\item{X}{（必须）数字向量。样地林木横坐标。单位：m。}

\item{Y}{（必须）数字向量。样地林木纵坐标。单位：m。}

\item{Origin}{数字向量，长度2。样地原点的坐标c(X0, Y0)。默认为c(0, 0)。单位：m。}

\item{Range_xy}{（必须）数值向量，长度2，均大于0。矩形样地横坐标X和纵坐标Y方向的长度。单位m。}
}
\value{
1）输出结果为一个数据框。

包括输入的Data数据框的所有列，以及函数产生的新列。

2）新输出变量名：

newX：调整后的新X坐标。

newY：调整后的新Y坐标。
}
\description{
Coord_Move()通过给予重复坐标一个随机的抖动来解决坐标重复问题。
}
\details{
1）Coord_Move()会产生新列newX和newY。

2）Coord_Move()中得到的新坐标与重复坐标的距离在0.1~0.3m之间，相对于重复坐标的方位是随机的。

3）当各样地均不存在重复坐标时，会输出“There are no duplicated coordinates!”的提醒，此时将不会再输出新数据框。

4）可对多个样地组成的数据框进行重复标签批量处理。
}
\examples{
## 加载内置数据
data(ForestStatTool)
##rawdata是一个存在重复坐标的样地原始数据框
##Coord_Move()通过给予重复坐标一个随机的抖动来解决坐标重复问题。
newdata1 <- Coord_Move(Data = rawdata, Plot = rawdata$Plot, X=rawdata$X, Y=rawdata$Y, Range_xy = c(10,10))
newdata1
##检查newdata1中的newX和newY中是否还存在重复坐标
Coord_Move(Data = newdata1, Plot = newdata1$Plot, X=newdata1$newX, Y=newdata1$newY, Range_xy = c(10,10))
#输出“There are no duplicated coordinates!”提醒，说明所有样地中的样木都不存在重复坐标
}
\author{
Chaofan Zhou: cfzhou2021@163.com

Xiao He: hexiao@ifrit.ac.cn

Guangshuang Duan: oliverdgs@163.com

Jingning Shi: shijingning@hebau.edu.cn

Xiangdong Lei: xdlei@ifrit.ac.cn
}
