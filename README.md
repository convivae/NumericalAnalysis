# NumericalAnalysis
数值分析（第四版） 颜庆津 书中部分算法的代码实现

## 第二章 线性方程组的解法

> input: “../2.txt”，增广矩阵

- 直接法
  - 高斯消去法
    - [ ] ~~顺序高斯消去法~~
    - [x] 高斯若当消去法
    - [x] 列主元素 Gauss 消去法
  - 直接三角分解法（A = LU）
    * [x] Doolittle 分解法
    * [x] Crout 分解法
    * [x] 选主元的 Doolittle 分解法
    * [ ] ~~三对角矩阵的追赶法~~
    * [ ] ~~平方根法~~
  - 矩阵的条件数
- 迭代法

  - [x] Jacobi 迭代法
  - [x] Gauss Seidel 迭代法
  - [x] SOR 迭代法

## 第三章 矩阵的特征值与特征向量的计算

> input: “../3.txt”，方阵

- 幂法及反幂法
  - [x] 幂法

    > 用于计算矩阵主特征值 (矩阵按模最大的特征值)及对应特征向量的迭代方法

  - [x] 反幂法

    - [x] 计算矩阵按模最小的特征值与其对应的特征向量
    - [x] 求一个给定近似特征值对应的特征向量

- Jacobi 方法

  - [ ] ~~平面旋转矩阵~~

  - [x] Jacobi 算法

    > 只适用于实对称方阵，求出所有特征值和特征向量

- QR 算法

  - [x] Householder 变换

    > 求解初等反射矩阵( Householder 矩阵)，使得将任意矩阵正交相似约化为对称三角阵

  - [x] QR 分解
  
    > 利用 Householder 矩阵（镜面映射矩阵）对矩阵 A 做 QR 分解，A = QR
  >
    > 其中，Q是正交矩阵，R是上三角矩阵

  - [x] QR 方法求矩阵的特征值（对应的特征向量可利用反幂法得到）
  
    - [x] QR 分解部分利用 Householder 矩阵直接得到
    - [x] QR 分解部分利用平面旋转矩阵得到
    - [x] 带原点位移的 QR 分解
    - [x] 带双步位移的 QR 分解

## 第四章 非线性方程与非线性方程组的迭代解法

> input: 迭代函数的函数指针

- 非线性方程的迭代解法

    * [x] 方程求根与二分法
    * [x] 简单迭代法
      - 简单迭代法的存在性与收敛性
      - 局部收敛性与收敛阶
    * [x] 迭代收敛的加速方法
      * [x] Steffensen 迭代法
    * [x] Newton 法
      - 单根的牛顿迭代法
        - [x] Newton 下山法
        - [ ] ~~简化牛顿法~~
    * [x] 割线法
      * [x] 割线法及其收敛性
      * [x] 单点割线法

- 非线性方程组的迭代解法

    * [x] 简单迭代法
    * [x] Newton 法

## 第五章 插值与逼近

> input: “../5.txt”，共两行， X<sub>i</sub> 和 Y<sub>i</sub>

* [x] Lagrange 插值
* [x] Newton 插值
* [ ] Hermite 插值
* [x] 样条插值
  * [x] 三弯矩法求三次样条插值
* [ ] 三角插值与快速 Fourier 变换
* [x] 正交多项式
  * [x] Legendre 多项式
  * [x] Chebyshev 多项式
  * [x] Laguerre 多项式
  * [x] Hermite 多项式
* [ ] 最佳平方逼近
* [ ] 曲线拟合与曲面拟合
  * [ ] 最小二乘法求拟合曲线
  * [ ] 最小二乘法求拟合曲面

## 第六章 数值积分

> input: 待积函数的函数指针

* [ ] Newton-Cotes求积公式
  * [ ] 梯形公式
  * [ ] Simpson 公式
  * [ ] Simpson 3/8 公式
  * [ ] Cotes 公式
* [ ] 复化求积法
  * [ ] 复化梯形公式
  * [ ] 复化 Simpson 公式
* [x] Romberg 积分法
* [x] Gauss 型求积公式
  * [x] Gauss-Legendre 求积公式
  * [x] Gauss-Laguerre 求积公式
  * [x] Gauss-Hermite 求积公式
  * [x] Gauss-Chebyshev 求积公式