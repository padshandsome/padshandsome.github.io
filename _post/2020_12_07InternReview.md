### AI and DA Intern in Shenzhen(Shutang Co)
## 一、 基础生物信息知识
[PPG](https://blog.csdn.net/wxq_1993/article/details/104022784)

采集途径：手腕、手指（光传感器）

PPG特征提取：
1. RR间期(RR interval)
峰值检测(peak detection) -> np.diff(peak) 

2. [HRV metrics](https://github.com/Aura-healthcare/hrvanalysis/blob/master/LICENSE)

根据Lancet一篇文章进行的PPG处理：

[PPG处理](_post/1607325480(1).png)


信号特点：易受动作影响，处理时要特别去掉motion artifact。

ECG

采集途径：手腕、手指

## 二、信号处理技术：
# 1. 滤波器
[Savitzky-Golay filter](https://blog.csdn.net/qq_20823641/article/details/51537461) 窗口越大越平滑。
普通滤波器：阻带频率计算=预期阻带频率/奈奎斯特频率=预期阻带频率/(采样频率/2)


## 三、机器学习
1. Adaboost

2. Random Forest


