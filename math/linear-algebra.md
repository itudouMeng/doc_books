# 向量空间

## 子空间(subspace)

矩阵的四个基本子空间：行空间(row space)、列空间(column space)、零空间(null space)、左零空间(left null space)。

行空间：矩阵行向量组成的线性空间，记为$C(A^T)$​​；

列空间：矩阵列向量组成的线性空间，记为$C(A)$；

零空间：表达式$Ax=0$​的解向量组成的线性空间，记为$N(A)$​；

左零空间：表达式$A^Tx=0$​​的解向量组成的线性空间，记为$N(A^T)$​​​。



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



矩阵行空间的维度和零空间的维度组成，即有表达式：
$$
\begin{aligned}
rank(A)=rank(A^T)=r\\
dim(C(A^T))+dim(N(A))=n\\
dim(C(A))+dim(N(A^T))=m
\end{aligned}
$$

<iframe src=https://www.bilibili.com/video/BV1bb411H7JN?p=16
        name="超赞的线性代数讲义"
        class="iframe"
        width=100% 
        height=600px
        scrolling="auto"
        sandbox="allow-scripts allow-forms allow-same-origin">
</iframe>


## 旋转矩阵

## 投影矩阵





## 反射矩阵

取直线$l$​​上单位方向向量$u=[\cos\theta,\sin\theta]^T$​​和单位法向向量$v=[\sin\theta,-\cos\theta]^T$​​作为一组基，那么这组基关于直线$l$对称后的基为$u'=u=[\cos\theta,\sin\theta]^T$​​和$v'=v=[-\sin\theta,\cos\theta]^T$​​。​​就有表达式：
$$
\begin{aligned}
P\left[\begin{matrix}
u&v
\end{matrix}\right]=
\left[\begin{matrix}
u'&v'
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
\left[\begin{matrix}
x'\\y'
\end{matrix}\right]
=\left[\begin{matrix}
\cos2\theta&\sin2\theta\\
\sin2\theta&-\cos2\theta
\end{matrix}\right]
\left[\begin{matrix}
x\\y
\end{matrix}\right]
$$
$P=\left[\begin{matrix}\cos2\theta&\sin2\theta\\\sin2\theta&-\cos2\theta\end{matrix}\right]$​​​​​​也是一个正定矩阵​，满足：(1)$P^{-1}=P^{T}$​，​(2)$P^2=P$​，(3)$P^TP=I$​



![refelction-matrix](linear-algebra.assets/refelction-matrix.svg)

