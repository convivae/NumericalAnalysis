# NumericalAnalysis
数值分析（第四版） 颜庆津 书中部分算法的代码实现

## 第二章 线性方程组的解法

* [x] 列主元素 Gauss 消去法
* [x] 高斯若当消去法
* [x] 三角分解法（A = LU）
  * [x] Doolittle 分解法
  * [x] Crout 分解法
  * [x] 选主元的 Doolittle 分解法
* [x] 迭代法
  * [x] Jacobi 迭代法
  * [x] Gauss Seidel 迭代法
  * [x] SOR 迭代法

## 第三章 矩阵的特征值与特征向量的计算

* [x] 幂法
  
  > 用于计算矩阵主特征值 (矩阵按模最大的特征值)及对应特征向量的迭代方法
  
* [x] 反幂法

  * [x] 计算矩阵按模最小的特征值与其对应的特征向量
  * [x] 求一个给定近似特征值对应的特征向量

* [x] Jacobi 方法

  > 只适用于实对称方阵，求出所有特征值和特征向量

* [x] QR 分解

  > 利用 Householder 矩阵（镜面映射矩阵）对矩阵 A 做 QR 分解，A = QR
  >
  > 其中，Q是正交矩阵，R是上三角矩阵

* [x] QR 方法求矩阵的特征值（对应的特征向量可利用反幂法得到）

  * [x] QR 分解部分利用 Householder 矩阵直接得到
  * [x] QR 分解部分利用平面旋转矩阵对矩阵得到
  * [x] 带原点位移的 QR 分解
  * [x] 带双步位移的 QR 分解

## 第四章 非线性方程与非线性方程组的迭代解法

### 4.1 非线性方程的迭代解法

* [ ] 简单迭代法
* [ ] Steffensen 迭代法
* [ ] Newton 法
* [ ] 割线法

### 4.2 非线性方程组的迭代解法

* [ ] 简单迭代法
* [ ] Newton 法