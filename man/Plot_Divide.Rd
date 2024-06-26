% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plot_Divide.R
\name{Plot_Divide}
\alias{Plot_Divide}
\title{样地分割}
\usage{
Plot_Divide(Data, Plot, X, Y, Num_xy = NULL, Length_xy = NULL, Origin = c(0,0), Range_xy)
}
\arguments{
\item{Data}{（必须）数据框。需要分割的样地数据框。}

\item{Plot}{（必须）向量，建议为字符型。样地号。}

\item{X}{（必须）数值向量。样地林木横坐标。单位：m。}

\item{Y}{（必须）数值向量。样地林木纵坐标。单位：m。}

\item{Num_xy}{（可选）数值向量，长度2。均大于1。横纵坐标方向划分小样地的个数。}

\item{Length_xy}{（可选）数值向量，长度2。均大于0且小于横纵坐标长度的一半。横纵坐标方向划分小样地的长度。单位：m。}

\item{Origin}{数值向量，长度2。样地原点的坐标c(X0, Y0)。默认为c(0, 0)。}

\item{Range_xy}{（必须）数值向量，长度2。均大于0。矩形样地横坐标X和纵坐标Y方向的长度。单位：m。}
}
\value{
1）输出结果为一个数据框。

包括输入的Data数据框的所有列，以及函数产生的新列。

2）新输出变量名：

subplot：分割后的新样地名，格式为“原样地名-小样地编号”。

subX：小样地的新X坐标，范围（0，Length_X）。单位：m。

subY：小样地的新X坐标，范围（0，Length_Y）。单位：m。
}
\description{
此函数实现了将矩形样地划分为面积相等的若干小样地的功能。可按照坐标X和Y方向上设置的数量或长度划分小样地。
}
\details{
1）可对多个样地组成的数据框进行样地分割的批量处理。

2）参数Length_X和Length_Y决定划分的小样地的大小，参数Num_X和Num_Y决定划分的小样地的数量。两组参数可以根据Length_X×Num_X = Xrange和Length_Y×Num_Y = Yrange相互推导，可以两组同时输入，或至少输入一组。当两组同时输入时，Length_X×Num_X和Length_Y×Num_Y分别需要小于或等于Xrange和Yrange，否则，将以参数Num_X和Num_Y为主，推导出Length_X和Length_Y。

3）小样地的编号规则：

当横坐标和纵坐标方向划分的小样地个数均小于等于10个时，采用列号（0-9）+行号（0-9）的2位数小样地编号。如34表示从左下角开始的第4列第5行的小样地。

当横坐标或纵坐标方向划分的小样地个数超过10个时，采用列号（00-99）+行号（00-99）的4位数小样地编号。如0706表示从左下角开始的第8列第7行的小样地。
}
\examples{
## 加载内置数据
data(ForestStatTool)
#bigplot是一个100m×100m的大样地数据框，用于分割小样地
#横纵坐标各划分5份，结果为划分后的25个20m×20m的小样地的合并数据框。
#输入方式
subplot1 <- Plot_Divide(Data = bigplot, Plot = bigplot$Plot, X = bigplot$X, Y = bigplot$Y,
                        Num_xy = c(5,5), Range_xy = c(100,100))
head(subplot1)
levels(factor(subplot1$subplot))
#一样的结果
subplot2 <- Plot_Divide(Data = bigplot, Plot = bigplot$Plot, X = bigplot$X, Y = bigplot$Y,
                        Length_xy = c(20,20), Range_xy = c(100,100))
head(subplot2)
levels(factor(subplot2$subplot))
#横坐标划分5份纵坐标各划分4份，结果为划分后的16个15m×20m的小样地的合并数据框
subplot3 <- Plot_Divide(Data = bigplot, Plot = bigplot$Plot, X = bigplot$X, Y = bigplot$Y,
                        Num_xy = c(5,4), Length_xy = c(15,20), Range_xy = c(100,100))
head(subplot3)
levels(factor(subplot3$subplot))
}
\author{
Chaofan Zhou: zhouchaofan@ifrit.ac.cn

Xiao He: hexiao@ifrit.ac.cn

Guangshuang Duan: oliverdgs@163.com

Jingning Shi: shijingning@hebau.edu.cn

Xiangdong Lei: xdlei@ifrit.ac.cn
}
