# 矩阵乘法

# 向量空间

## 子空间

矩阵的四个基本子空间(subspace)：行空间(row space)、列空间(column space)、零空间(null space)、左零空间(left null space)。

- **行空间：矩阵行向量组成的线性空间，记为$C(A^T)$​​；**
- **列空间：矩阵列向量组成的线性空间，记为$C(A)$；**
- **零空间：表达式$Ax=0$​的解向量组成的线性空间，记为$N(A)$​；**
- **左零空间：表达式$A^Tx=0$​​​的解向量组成的线性空间，记为$N(A^T)$​​​​​​​。**

当$m\times n$​阶矩阵A的秩$rank(A)=r$​,可以通过高斯消元法得到：
$$
\begin{aligned}
A\vec x=\bold 0 &\rArr
\left[\begin{matrix}
I&F\\
\bold 0&\bold 0
\end{matrix}\right]
\vec x=\bold 0 

\\&\rArr \vec x=
\left[\begin{matrix}
-F\\I  
\end{matrix}\right]
\end{aligned}
$$

矩阵的行空间与零空间的维数之和为矩阵的列数$n$​​，并且两个子空间相互正交:
$$
dim(C(A^T))+dim(N(A))=n\\
C(A^T)\cdot N(A)=0
$$
矩阵的列空间与左零空间维数之和为矩阵的行数$m$，并且两个子空间相互正交：
$$
dim(C(A))+dim(N(A^T))=m\\
C(A)\cdot N(A^T)=0
$$
![4-subspaces](linear-algebra.assets/4-subspaces.png)


<iframe src=https://www.bilibili.com/video/BV1bb411H7JN?p=16
        name="超赞的线性代数讲义"
        class="iframe"
        width=100% 
        height=600px
        scrolling="auto"
        sandbox="allow-scripts allow-forms allow-same-origin">
</iframe>



## 投影矩阵

向量投影变换是指向量在矩阵列空间上的投影，记投影变换矩阵为$P$​​​​​​，即投影变换为:$\boldsymbol b\rightarrow P\boldsymbol b$​​​​​​。记向量$\boldsymbol b$​​​​​​在矩阵$A$​​​​​​的列空间上的投影为$\boldsymbol p$​​​​​​，那么二者的差向量$\boldsymbol b-\boldsymbol p$​​​​​​就应当垂直于投影$\boldsymbol p$​​​​​​​(投影的定义，这样保证投影向量是“最小的”)，所以就有表达式：
$$
\boldsymbol p\cdot(\boldsymbol b-\boldsymbol p)=\boldsymbol 0
\Rightarrow (P\boldsymbol b)^T\cdot(\boldsymbol b-P\boldsymbol b)=\boldsymbol 0
$$

显然从上面的表达式中无法求得$P$​​，再补充一条件，投影$\boldsymbol p$​​在矩阵的矩阵$A$​​的列空间上(投影$\boldsymbol p$​​可以由矩阵$A$​​的列向量线性组合而来)，即满足：
$$
\boldsymbol p=A\hat{\boldsymbol x}
$$
代入式中，然后我们就能够得到：
$$
A^T(\boldsymbol b-A\hat{\boldsymbol x})=\boldsymbol 0
\Rightarrow A^T A\hat{\boldsymbol x}=A^T\boldsymbol b
$$


![projection-onto-column-space](linear-algebra.assets/projection-onto-column-space.png)

实际上我们可以这么来理解：差向量$\boldsymbol e=\boldsymbol b-\boldsymbol p=\boldsymbol b-A\hat{\boldsymbol x}$​垂直于矩阵$A$​的列空间，就可以直接得到上式。当矩阵$A$​的列线性无关时，对上面的表达式化简，我们就可以得到投影变换为：
$$
\boldsymbol p=A(A^T A)^{-1}A^T\boldsymbol b
$$
其中，矩阵$P=A(A^T A)^{-1}A^T$​就是投影矩阵。显然，当矩阵的列空间是一条直线时，也是满足投影变换的，这时候我们就有：
$$
\boldsymbol p=\boldsymbol a\hat x=\frac{\boldsymbol a\boldsymbol a^T}{\boldsymbol a^T\boldsymbol a}\boldsymbol b
$$

![projection-onto-line](linear-algebra.assets/projection-onto-line.png)

实际应用中，利用矩阵乘法和向量数乘性质，能够较快地得到向量投影计算式：
$$
\boldsymbol p=\frac{\boldsymbol a\boldsymbol a^T}{\boldsymbol a^T\boldsymbol a}\boldsymbol b=\boldsymbol a\frac{\boldsymbol a^T\boldsymbol b}{\boldsymbol a^T\boldsymbol a}=\hat x\boldsymbol a
$$

## 旋转矩阵

记平面向量关于坐标原点旋转$\theta$​角后的旋转变换为：$\boldsymbol b\rightarrow Q\boldsymbol b$​​​。

![rotate-matrix](linear-algebra.assets/rotate-matrix.png)

取坐标系下一组单位基：$\boldsymbol u=[1,0]^T$​​和$\boldsymbol v=[0,1]^T$​​，很容易得到旋转之后在原坐标系下变为：$\boldsymbol u'=[\cos\theta,\sin\theta]^T$​​和$\boldsymbol v=[-\sin\theta,\cos\theta]^T$​。设旋转矩阵为$Q$，那么，就有：​​​​
$$
Q\left[\begin{matrix}
\boldsymbol u&\boldsymbol v
\end{matrix}\right]
=\left[\begin{matrix}
\boldsymbol u'&\boldsymbol v'
\end{matrix}\right]\Rightarrow
Q\left[\begin{matrix}
1&0\\
0&1
\end{matrix}\right]=
\left[\begin{matrix}
\cos\theta&-\sin\theta\\
\sin\theta&\cos\theta
\end{matrix}\right]
$$
即可求矩阵$Q$​得表达式为：


$$
Q=\left[\begin{matrix}
\cos\theta&-\sin\theta\\
\sin\theta&\cos\theta
\end{matrix}\right]
$$

相似地，在三维空间中我们可以得到，绕$x$轴旋转$\alpha$角的旋转矩阵为：
$$
Q_{\alpha}=\left[\begin{matrix}
1&0&0\\
0&\cos\alpha&-\sin\alpha\\
0&\sin\alpha&\cos\alpha
\end{matrix}\right]
$$
绕$y$​轴旋转$\beta$​角的旋转矩阵为：
$$
Q_{\beta}=\left[\begin{matrix}
\cos\beta&0&\sin\beta\\
0&1&0\\
-\sin\beta&0&\cos\beta
\end{matrix}\right]
$$
绕$z$​​轴旋转$\gamma$​​​角的旋转矩阵为：
$$
Q_{\gamma}=\left[\begin{matrix}
\cos\gamma&-\sin\gamma&0\\
\sin\gamma&\cos\gamma&0\\
0&0&1
\end{matrix}\right]
$$
那么三维空间中旋转矩阵可以表示为：
$$
Q=Q_{\alpha}Q_{\beta}Q_{\gamma}
$$

## 反射矩阵

记向量关于直线$l$​​​​​​反射变换矩阵为$P$​​​​​​，即反射变换为：$\boldsymbol b\rightarrow P\boldsymbol b$​​​​​​。取直线$l$​​​​​​上单位方向向量$\boldsymbol u=[\cos\theta,\sin\theta]^T$​​​​​​和垂直于直线的单位法向向量$\boldsymbol v=[\sin\theta,-\cos\theta]^T$​​​​​​作为一组基，关于直线$l$​​​​​​对称后的基为$\boldsymbol u'=\boldsymbol u=[\cos\theta,\sin\theta]^T$​​​​​​和$\boldsymbol v'=-\boldsymbol v=[-\sin\theta,\cos\theta]^T$​​​​​​​。就有表达式：
$$
\begin{aligned}
P\left[\begin{matrix}
\boldsymbol u&\boldsymbol v
\end{matrix}\right]=
\left[\begin{matrix}
\boldsymbol u'&\boldsymbol v'
\end{matrix}\right]&\Rightarrow
P\left[\begin{matrix}
\cos\theta&\sin\theta\\
\sin\theta&-\cos\theta\\
\end{matrix}\right]=
\left[\begin{matrix}
\cos\theta&-\sin\theta\\
\sin\theta&\cos\theta\\
\end{matrix}\right]
\\&\Rightarrow
P=\left[\begin{matrix}\cos2\theta&\sin2\theta\\\sin2\theta&-\cos2\theta\end{matrix}\right]
\end{aligned}
$$
那么，平面坐标系下关于直线对称的矩阵变换是：
$$
\begin{bmatrix}
x'\\y'
\end{bmatrix}
=\begin{bmatrix}
\cos2\theta&\sin2\theta\\
\sin2\theta&-\cos2\theta
\end{bmatrix}
\begin{bmatrix}
x\\y
\end{bmatrix}
$$
$P=\begin{bmatrix}\cos2\theta&\sin2\theta\\\sin2\theta&-\cos2\theta\end{bmatrix}$​​​​​​​​也是一个正交矩阵​，满足：(1)$P^{-1}=P^{T}$​​​​​，​(2)$P^2=I$​

## 格拉姆-施密特正交化

格拉姆-施密特正交化(Gram-Schmidt)方法一般用于一组标准规范基的构造，原理主要是利用原向量与投影向量的差向量正交于投影向量的性质。设有$n$​​​个线性无关的向量：$\boldsymbol a_{1},\boldsymbol a_{2},\cdots,\boldsymbol a_{n}$​​​，通过正交化方法最后形成一组标准规范基：$\boldsymbol q_{1},\boldsymbol q_{2},\cdots,\boldsymbol q_{n}$​​​。其正交化过程可以表述为：

- step1:选取第一个元素，进行标准化​

$$
\boldsymbol q_{1}=\frac{\boldsymbol a_{1}}{\Vert\boldsymbol a_{1}\Vert}
$$

- step2:对第二个元素进行正交化，并且标准化
$$
\begin{cases}
\boldsymbol q_{2}=\boldsymbol a_{2}-\frac{\boldsymbol q_{1}\boldsymbol q_{1}^T}{\boldsymbol q_{1}^T\boldsymbol q_{1}}\boldsymbol a_{2}\\
\boldsymbol q_{2}=\frac{\boldsymbol q_{2}}{\Vert\boldsymbol q_{2}\Vert}
\end{cases}
$$

- step3:利用投影矩阵的性质，对第三个元素进行正交化和标准化

$$
\begin{cases}
\boldsymbol q_{3}=\boldsymbol a_{3}-\frac{\boldsymbol q_{1}\boldsymbol q_{1}^T}{\boldsymbol q_{1}^T\boldsymbol q_{1}}\boldsymbol a_{3}-\frac{\boldsymbol q_{2}\boldsymbol q_{2}^T}{\boldsymbol q_{2}^T\boldsymbol q_{2}}\boldsymbol a_{3}\\
\boldsymbol q_{3}=\frac{\boldsymbol q_{3}}{\Vert\boldsymbol q_{3}\Vert}
\end{cases}
$$

- ...
- step n:依次类推，可以得到第$n$​​个规范基的表达式：

$$
\begin{cases}
\boldsymbol q_{n}=\boldsymbol a_{n}-\frac{\boldsymbol q_{1}\boldsymbol q_{1}^T}{\boldsymbol q_{1}^T\boldsymbol q_{1}}\boldsymbol a_{n}-\frac{\boldsymbol q_{2}\boldsymbol q_{2}^T}{\boldsymbol q_{2}^T\boldsymbol q_{2}}\boldsymbol a_{n}-\cdots-\frac{\boldsymbol q_{n-1}\boldsymbol q_{n-1}^T}{\boldsymbol q_{n-1}^T\boldsymbol q_{n-1}}\boldsymbol a_{n}\\
\boldsymbol q_{n}=\frac{\boldsymbol q_{n}}{\Vert\boldsymbol q_{n}\Vert}
\end{cases}
$$

将正交化结果用矩阵形式表示，就有：
$$
A=QQ^TA=Q(Q^TA)=QR
$$

即有下式：
$$
\begin{bmatrix}
\\
\boldsymbol a_1&\boldsymbol a_2&\cdots&\boldsymbol a_n\\
\\
\end{bmatrix}=
\begin{bmatrix}
\\
\boldsymbol q_1&\boldsymbol q_2&\cdots&\boldsymbol q_n\\
\\
\end{bmatrix}
\begin{bmatrix}
\boldsymbol q_{1}^T\boldsymbol a_{1}&\boldsymbol q_{1}^T\boldsymbol a_{2}&\cdots&\boldsymbol q_{1}^T\boldsymbol a_{n}\\
&\boldsymbol q_{2}^T\boldsymbol a_{2}&\cdots&\boldsymbol q_{2}^T\boldsymbol a_{n}\\
&&\ddots&\vdots\\
&&&\boldsymbol q_{n}^{T}\boldsymbol a_{n}
\end{bmatrix}
$$

上面的表达式就是我们通常意义上所说的QR分解，即$A=QR$​。​

## 最小二乘法

最小二乘法(least squares approximation)在统计学和工程计算中占有重要的地位，通过最小化误差来进行数据拟合。实际中，输入数据本身是存在误差的，精确解通常是无法得到的。一般情况下，对多组输入数据，我们可以得到一个矩阵方程：
$$
Ax=\boldsymbol b
$$
当方程数量大于未知数数量时（系数矩阵$A$的行数大于列数：$m>n$），原方程不一定存在解，这个时候最小二乘法就能求得一个近似解$\hat x$，某种程度上其是“最优解“，几何意义上平方差$error(x-\hat x)$​​最小。

设最优解是$\hat x$​,那么最优解$\hat x$​在系数矩阵$A$​的列空间上的投影就为$\boldsymbol p=A\hat x$​，向量$\boldsymbol b$与投影$\boldsymbol p$的差向量就为$\boldsymbol e=\boldsymbol b-\boldsymbol p=\boldsymbol b-A\hat x$​​,差向量应当垂直于系数矩阵$A$​的列空间。我们就可以得到：
$$
A^T\boldsymbol e=0\rightarrow A^TA\hat x=A^T\boldsymbol b
$$

# 矩阵特征值和特征向量

设$n$阶矩阵$A$有$n$个特征值(eigenvalues)，记为：$\lambda_{1},\lambda_{2},\cdots,\lambda_{n}$​，有：
$$
det(A-\lambda I)=\begin{vmatrix}
a_{11}-\lambda&&&\\
&a_{22}-\lambda&&\\
&&\ddots&\\
&&&a_{nn}-\lambda
\end{vmatrix}
=(\lambda_{1}-\lambda)(\lambda_{2}-\lambda)\cdots(\lambda_{n}-\lambda)=0
$$
根据方程根与系数之间的关系 ，展开可得到矩阵特征值与矩阵元素之间的关系，其中比较重要的有下面两条。

- 矩阵的特征值之和等于该矩阵的迹(trace)，即矩阵的对角元素之和：

$$
\sum_{i=1}^{n} \lambda_{i}=trac(A)=\sum_{i=1}^{n} a_{ii}
$$

- 矩阵的特征值之积等于该矩阵的行列式(determinate)值，有：

$$
\prod_{i=1}^{n}\lambda_{i}=det(A)
$$

## 矩阵对角化

矩阵对角化(diagnalize)在处理矩阵的幂乘(power of matrix)时非常有用，当矩阵有n个线性无关特征向量时，矩阵可对角化。

根据矩阵特征值的定义，有：
$$
A\boldsymbol x_{1}=\lambda_{1}\boldsymbol x_{1}\\
A\boldsymbol x_{2}=\lambda_{2}\boldsymbol x_{2}\\
\cdots\\
A\boldsymbol x_{n}=\lambda_{n}\boldsymbol x_{n}
$$
将上述结果写成矩阵形式，我们就可以得到：
$$
AX=X\Lambda
$$
其中，特征向量矩阵$X=[\boldsymbol x_{1},\boldsymbol x_{2},\cdots,\boldsymbol x_{n}]$​，对角矩阵$\Lambda=diag(\lambda_{1},\lambda_{2},\cdots,\lambda_{n})$​。​

因为特征向量线性无关，即矩阵$X$可逆，进一步，我们可以得到矩阵的对角化结果：
$$
X^{-1}AX=\Lambda\Leftrightarrow A=X\Lambda X^{-1}
$$
矩阵的幂乘就可以简单地表示为：
$$
A^{k}=\begin{matrix}
\underbrace{(X\Lambda X^{-1})(X\Lambda X^{-1})\cdots(X\Lambda X^{-1})}\\
n
\end{matrix}
=X\Lambda^{k}X^{-1}
$$

## 矩阵对角化判别规则

几何重数(GM=geometric multiplicity)指对应特征值$\lambda$的线性无关的特征向量的数量，几何重数是$A-\lambda I$​的零空间的维度，即：
$$
GM=dim(N(A-\lambda I))
$$
代数重数(AM=algebric multiplicity)指特征多项式$det(A-\lambda I)=0$​​根的数量，一般地，代数重数是矩阵的阶数，即：
$$
AM=n
$$
通常有关系式：$GM\leq AM=n$​​。当几何重数小于代数重数时，矩阵不可对角化；当几何重数等代数重数时，矩阵可对角化。

## 凯莱-哈密顿定理

凯莱-哈密顿定理(Caylay-Hamilton Theorem)描述了矩阵特征值与矩阵的一种关系。记矩阵的特征多项式为$f(\lambda)=(\lambda-\lambda_1)(\lambda-\lambda_2)\cdots(\lambda-\lambda_n)$​，那么有：
$$
f(A)=(A-\lambda_1 I)(A-\lambda_2 I)\cdots(A-\lambda_n I)=\boldsymbol O
$$
当矩阵可对角化时，就有：
$$
\begin{aligned}
(A-\lambda_1 I)\cdots(A-\lambda_n I)&=(X\Lambda X^{-1}-X\lambda_1X^{-1})\cdots(X\Lambda X^{-1}-X\lambda_nX^{-1})
\\&=X(\Lambda-\lambda_1 I)X^{-1}\cdots X(\Lambda-\lambda_n I)X^{-1}\\&=X\begin{bmatrix}\lambda_1-\lambda_1&&\\&\ddots&\\&&\lambda_n-\lambda_1\end{bmatrix}\cdots
\begin{bmatrix}\lambda_1-\lambda_n&&\\&\ddots&\\&&\lambda_n-\lambda_n\end{bmatrix}X^{-1}
\\&=XOX^{-1}=O
\end{aligned}
$$
当矩阵不可对角化时，亦有关系式成立。

## 斐波那契数列

对于递推关系式$\boldsymbol u_{k+1}=A\boldsymbol u_{k}$和初始值$\boldsymbol u_{0}$​，由递推关系我们可以得到：
$$
\boldsymbol u_{k}=A^{k}\boldsymbol u_{0}
$$
当矩阵$A$可对角化时，我们对上式化简就可以得到：
$$
\boldsymbol u_{k}=X\Lambda^{k}X^{-1}\boldsymbol u_{0}\overset{let\ X^{-1}\boldsymbol u_{0}=\boldsymbol c}{=}c_{1}\lambda_{1}^{k}\boldsymbol x_{1}+c_{2}\lambda_{2}^{k}\boldsymbol x_{2}+\cdots+c_{n}\lambda_{n}^{k}\boldsymbol x_{n}
$$
上述的关键就是求得系数矩阵$A$的特征值和特征向量。

我们知道，斐波那契数列(Fibonacci Numbers)可以表示为：
$$
F_{n+2}=F_{n+1}+F_{n}\\
F_{1}=0,F_{2}=1
$$
那么，对于斐波那契数列，我们就有一般的求解过程：

- step1：对递推关系式构造系数矩阵

$$
\begin{aligned}
F_{n+2}&=F_{n+1}+F_{n}\\
F_{n+1}&=F_{n+1}
\end{aligned}
$$

记$\boldsymbol u_{k}=\begin{bmatrix}F_{n+1}\\F_{n}\end{bmatrix}$​，​​那么上式用矩阵的形式就可以表式为：
$$
\begin{aligned}
\boldsymbol u_{k+1}&=\begin{bmatrix}
1&1\\
1&0
\end{bmatrix}
\boldsymbol u_{k}\\
\boldsymbol u_{0}&=
\begin{bmatrix}
1\\0
\end{bmatrix}
\end{aligned}
$$

- step2：求系数矩阵的特征值和特征向量

$$
det(A-\lambda I)=\begin{vmatrix}
1-\lambda&0\\
1&-\lambda
\end{vmatrix}=\lambda^{2}-\lambda-1=0
$$

求得特征值为：$\lambda_{1}=\frac{1+\sqrt{5}}{2}$​​​和$\lambda_{2}=\frac{1-\sqrt{5}}{2}$​​​。特征向量为：$\boldsymbol x_{1}=\begin{bmatrix}\lambda_{1}\\1\end{bmatrix}$​​​和$\boldsymbol x_{2}=\begin{bmatrix}\lambda_{2}\\1\end{bmatrix}$​​​​​。​

- step3：代入递推关系式，求得表达式

$$
\begin{aligned}
\boldsymbol u_{k}&=c_{1}\lambda_{1}^{k}\boldsymbol x_{1}+c_{2}\lambda_{2}^{k}\boldsymbol x_{2}
\\&=c_{1}\lambda_{1}^{k}\begin{bmatrix}\lambda_{1}\\1\end{bmatrix}+c_{2}\lambda_{2}^{k}\begin{bmatrix}\lambda_{2}\\1\end{bmatrix}
\\&=\begin{bmatrix}
c_{1}\lambda^{k+1}+c_{2}\lambda^{k+1}\\
c_{1}\lambda^{k}+c_{1}\lambda^{k}
\end{bmatrix}
\end{aligned}
$$

根据初始值$\boldsymbol u_{0}=\begin{bmatrix}1\\0\end{bmatrix}$​​​​，求得$c_{1}=\frac{1}{\lambda_{1}-\lambda_{2}}=\frac{\sqrt{5}}{5}$​​​​和$c_{2}=\frac{-1}{\lambda_{1}-\lambda_{2}}=-\frac{\sqrt{5}}{5}$​​​​​​。于是，我们就得到斐波那契数列的通项为：
$$
F_{n}=\frac{\sqrt{5}}{5}(\frac{1+\sqrt{5}}{2})^{n}-\frac{\sqrt{5}}{5}(\frac{1-\sqrt{5}}{2})^{n}
$$
-----------------------------------------***----------------------------------------------

当然了，我们也可以从矩阵的维度出发，得到：
$$
\boldsymbol u_{k}=X\Lambda^{k}X^{-1}\boldsymbol u_{0}\overset{let\ X^{-1}\boldsymbol u_{0}=\boldsymbol c}{=}X\Lambda\boldsymbol c
$$
于是就有：
$$
X^{-1}\boldsymbol u_{0}=\boldsymbol c\Rightarrow \boldsymbol u_{0}=X\boldsymbol c
$$
可以理解$\boldsymbol c$​为就是初始值$\boldsymbol u_{0}$​在特征向量矩阵$X$​​上的投影。所以step3可以替换为：

- step3'：求得初始值向量$u_{0}$在特征向量$X$上的投影坐标

$$
\boldsymbol u_{0}=c_{1}\boldsymbol x_{1}+c_{2}\boldsymbol x_{2}\Rightarrow
\begin{bmatrix}
1\\0
\end{bmatrix}
=c_{1}\begin{bmatrix}
\lambda_{1}\\1
\end{bmatrix}+
c_{2}\begin{bmatrix}
\lambda_{2}\\1
\end{bmatrix}
\Rightarrow c_{1}=\frac{1}{\lambda_{1}-\lambda_{2}}&c_{2}=\frac{-1}{\lambda_{1}-\lambda_{2}}
$$

带入式$\boldsymbol u_{k}=c_{1}\lambda_{1}^{k}\boldsymbol x_{1}+c_{2}\lambda_{2}^{k}\boldsymbol x_{2}$​​，于是就能得到斐波那契数列的通项。

-----------------------------------------***----------------------------------------------

## 矩阵指数

类似指数函数的泰勒展开式，矩阵指数(matrix exponent)$e^{At}$​​有定义式：
$$
e^{At}=I+At+\frac{(At)^{2}}{2!}+\frac{(At)^{3}}{3!}+\cdots=\sum_{k=1}^{\infin}\frac{(At)^{k}}{k!}
$$

- 导函数

  我们对变量$t$​​求导，有：

$$
A+A^{2}t+\frac{A^{3}t^{2}}{2!}+\cdots=Ae^{At}
$$

​	我们就有一个优美和简洁的结论，矩阵指数$e^{At}$的导函数与指数函数一般形式保持一致。

- 特征值和特征向量
记矩阵$A$​特征值为$\lambda$，特征向量为$\boldsymbol x$​​​，对矩阵指数乘以特征向量$\boldsymbol x$​就得到：

$$
(I+At+\frac{(At)^{2}}{2!}+\cdots)\boldsymbol x=(1+\lambda t+\frac{(\lambda t)^{2}}{2!}+\cdots)\boldsymbol x=e^{\lambda t}\boldsymbol x
$$

​	我们就有矩阵指数$e^{At}$​的特征值为$e^{\lambda t}$​，特征向量同为矩阵$A$​的特征向量$\boldsymbol x$​​。

- 对角化

  利用表达式$A=X\Lambda X^{-1}$，对矩阵指数进行对角化，有：

$$
\begin{aligned}
e^{At}&=I+X\Lambda tX^{-1}+\frac{(X\Lambda t X^{-1})(X\Lambda t X^{-1})}{2!}+\cdots\\
&=X\left[I+\Lambda t+\frac{(\Lambda t)^{2}}{2!}+\cdots\right]X^{-1}\\
&=Xe^{\Lambda t}X^{-1}\\
&=X\begin{bmatrix}
e^{\lambda_1}&&\\
&\ddots&\\
&&e^{\lambda_n}
\end{bmatrix}X^{-1}
\end{aligned}
$$

- 稳态

  当$\bold u(t)=e^{At}\boldsymbol u_0\rightarrow 0$​​​​时，系统最终趋于稳定状态。就有当矩阵$A$​的所有特征值实部$Re(\lambda)$小于$0$​时，满足条件$e^{At}\rightarrow 0$​。

## 常系数微分方程

记$\boldsymbol u$​​​​是时间$t$​​​​的函数：$\boldsymbol u=\boldsymbol u(t)$​​​​,初始值为$\boldsymbol u_{0}$。那么，对于常系数微分方程
$$
\begin{cases}
\begin{aligned}
\frac{d\boldsymbol u}{dt}&=A\boldsymbol u\\
\boldsymbol u_0&=\boldsymbol u(0)
\end{aligned}
\end{cases}
$$
由矩阵指数的定义，其矩阵形式的通解就为：
$$
\boldsymbol u=e^{At}\boldsymbol u_{0}
$$

1. 当矩阵$A$​​​​​可对角化时，我们就有:

$$
\begin{aligned}
   \boldsymbol u=e^{At}\boldsymbol u_{0}&=Xe^{\Lambda t}X^{-1}\boldsymbol u_{0}
   \\&\overset{let\ X^{-1}\boldsymbol u_{0}=\boldsymbol c}{=}Xe^{\Lambda t}\boldsymbol c
   \\&=c_{1}e^{\lambda_{1}t}\boldsymbol x_{1}+c_{2}e^{\lambda_{2}t}\boldsymbol x_{2}+\cdots
   \end{aligned}
$$

2. 当矩阵$A$​​​不可对角化时，表示矩阵的几何重数小于代数重数，即矩阵不具有$n$​​个线性无关的特征向量,存在重复的特征向量。记矩阵重复特征向量对应的特征值为$\lambda_{0}$​​，我们就有:

$$
\begin{aligned}
\boldsymbol u=e^{At}\boldsymbol u_{0}&=e^{(A-\lambda_{0}I)t+\lambda_{0}It}\boldsymbol u_{0}
\\&=e^{\lambda_{0}t}e^{A-\lambda_{0}I}\boldsymbol u_{0}
\\&=e^{\lambda_{0}t}\left[I+(A-\lambda_{0}I)t+\frac{[(A-\lambda_{0}I)t]^2}{2!}+\cdots\right]\boldsymbol u_{0}
\end{aligned}
$$

​	可以通过矩阵指数展开式进一步分析解的形式。

类似的，我们可以得到常系数微分方程的一般求解过程。

### 动力系统耦合方程

以下面的耦合方程来演示求解过程：
$$
\begin{cases}
\begin{aligned}
\frac{d\boldsymbol u_1}{dt}&=\boldsymbol u_1+\boldsymbol u_2\\
\frac{d\boldsymbol u_2}{dt}&=\boldsymbol u_1-\boldsymbol u_2
\end{aligned}
\end{cases}
$$
初始值为：$\boldsymbol u(0)=\begin{bmatrix}\boldsymbol u_1(0)\\\boldsymbol u_2(0)\end{bmatrix}=\begin{bmatrix}1\\0\end{bmatrix}$。

- step1:构造微分方程的矩阵表达形式

$$
\frac{d\boldsymbol u}{dt}=\begin{bmatrix}
1&1\\
1&-1
\end{bmatrix}\begin{bmatrix}
\boldsymbol u_{1}\\
\boldsymbol u_{2}
\end{bmatrix}=A\boldsymbol u
$$

- step2:求得系数矩阵的特征值和特征向量

  根据特征多项式$det(A-\lambda I)=0$​​​，求得特征值为：$\lambda_1=\sqrt2$​​​和$\lambda_2=-\sqrt2$​​​；特征向量为：$\boldsymbol x_1=\begin{bmatrix}1+\lambda_1\\1\end{bmatrix}$​​​和$\boldsymbol x_2=\begin{bmatrix}1+\lambda_2\\1\end{bmatrix}$​​​。

- step3:求得初始值在特征向量上的投影

$$
\begin{aligned}
\boldsymbol u_0=c_1\boldsymbol x_1+c_2\boldsymbol x_2
&\Rightarrow\begin{bmatrix}1\\0\end{bmatrix}=c_1\begin{bmatrix}1+\lambda_1\\0\end{bmatrix}+c_2\begin{bmatrix}1+\lambda_2\\0\end{bmatrix}
\\&\Rightarrow\begin{bmatrix}c_1\\c_2\end{bmatrix}=\begin{bmatrix}\frac{1}{\lambda_1-\lambda_2}\\\frac{1}{\lambda_2-\lambda_1}\end{bmatrix}=\begin{bmatrix}\frac{\sqrt2}{4}\\-\frac{\sqrt2}{4}\end{bmatrix}
\end{aligned}
$$

- step4:代入表达式，得到通解

$$
\begin{aligned}
\boldsymbol u=e^{At}\boldsymbol u_{0}&=c_{1}e^{\lambda_{1}t}\boldsymbol x_{1}+c_{2}e^{\lambda_{2}t}\boldsymbol x_{2}
\\&=\frac{\sqrt2}{4}e^{\sqrt2t}\begin{bmatrix}1+\sqrt2\\1\end{bmatrix}-\frac{\sqrt2}{4}e^{-\sqrt2t}\begin{bmatrix}1-\sqrt2\\1\end{bmatrix}
\\&=\begin{bmatrix}\frac{\sqrt2}{4}(1+\sqrt2)e^{\sqrt2t}-\frac{\sqrt2}{4}(1-\sqrt2)e^{-\sqrt2t}\\\frac{\sqrt2}{4}e^{\sqrt2t}-\frac{\sqrt2}{4}e^{-\sqrt2t}\end{bmatrix}
\end{aligned}
$$

### 二阶微分方程

对于二阶微分方程$m\ddot y+c\dot y+ky=0$​​​​​​，可以通过构造矩阵的形式来求解：
$$
\begin{cases}
\begin{aligned}
\ddot y&=-\frac{c}{m}\dot y-\frac{k}{m}y\\
\dot y&=\dot y
\end{aligned}
\end{cases}\Rightarrow
\begin{bmatrix}\ddot y\\\dot y\end{bmatrix}
=\begin{bmatrix}-\frac{c}{m}&-\frac{k}{m}\\1&0\end{bmatrix}
\begin{bmatrix}\dot y\\ y\end{bmatrix}
\Rightarrow \dot u=Au
$$

微分方程的详细解法参考：[differential equations and linear algebra]:https://math.mit.edu/dela

## 马尔可夫矩阵

马尔可夫矩阵在随机过程中有重要的应用。马尔可夫矩阵有定义：

1. 所有元素大于等于O
2. 各列向量元素之和为1

基本性质：

- 1是特征值
- 其余特征值绝对值小于1。当系统可以用马尔可夫矩阵描述时，其稳态是特征值1对应的那部分矢量：

$$
\begin{cases}
\begin{aligned}
\boldsymbol u&=\boldsymbol u(t)\\
\boldsymbol u_k&=A^k\boldsymbol u_0\to\lim_{k\to\infin}\boldsymbol u_k=c_1e^{t}\boldsymbol x_1
\end{aligned}
\end{cases}
$$

## <u>海森堡不确定性原则</u>

> $AB-BA=I$ can happen for infinite matricies with $A=A^T$ and $B=-B^T$​​.Then
> $$
> x^Tx=x^TABx-x^TBAx \leq  2\Vert Ax \Vert\Vert Bx \Vert
> $$
> Explain that last step by using the Schwarz inequality$|u^Tv|\leq \Vert u \Vert\Vert v \vert$.Then Heisenberg's inequality says that $\Vert Ax \Vert/\Vert x \Vert$ times $\Vert Bx \Vert/\Vert x \Vert$​is at least ½­ It is impossible to get the position error and momentum error both very small.



## 对称矩阵

