# 退化数据在可靠性中的应用

## 目的

在 R 语言环境中，关于随机过程在可靠性应用的案例与代码很少。这使得工程师在使用随机过程进行可靠性分析时，存在瓶颈。本项目是对一些经典的退化过程（维纳、逆高斯、伽马、泊松、指数扩散过程）进行学习并将其复现为 R 代码，最后制作成 R 包以及给出简单教程供大家学习和参考。


## 预计实现内容：

1. 学习各类经典退化过程，利用 bookdown 整理成一本笔记本；

2. 制作一个相关 R 包（Dagradation Process in Reliability, DPIR）；

3. 整理成可发表的论文，预计发表一篇 SCI 论文。


## 进度表

|序号|    进程     |    预计完成时间   |   真实完成时间    |    备注   |
|:----------:|:----------:|:----------:|:--------:|:--------:|
|1| 汇总退化过程在可靠性中应用案例与方法。||||
|2| 编写经典随机过程分析代码，包括（生成随机过程，统计推断，模型比较，可视化等）||||
|3| 制作成 R 包，并公开发布提供简单教程||||

## 备注

- 该项目与 [`退化数据的收集`] 和 [`退化过程综述`] 相结合，同步进行。需要一些合作者一起加入。

- 第二部分需要大量时间，还需要根据第一部分的学习进行细分。例如：

    -  具体做哪几个过程（维纳，伽马，逆高斯，指数扩散过程），参考综述《Remaining useful life estimation – A review on the statistical data driven approaches》？
    -  拓展模型用哪几个（随机效应，测量误差，多阶段，多失效机制），参考综述《Degradation data analysis and remaining useful life estimation: A review on Wiener-process-based methods
Remaining useful life estimation – A review on the statistical data driven approaches》？ 
    -  统计推断方法（极大似然，EM算法，贝叶斯方法以及非参数方法），是否需要使用拓展的方法？
    -  能否使用可以共用的统计推断方法？（MLE的话需要每个过程手写对数似然函数，工作量大！）
    -  新的推断方法是否可以引进（贝叶斯最新的方法，近似贝叶斯，合成似然？）
    -  案例分析（具体使用哪些典型的方法）







