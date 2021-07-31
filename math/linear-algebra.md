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

[超赞的线性代数讲义](https://www.bilibili.com/video/BV1bb411H7JN?p=16&spm_id_from=pageDriver ':include :type=iframe width=100% height=600px')

