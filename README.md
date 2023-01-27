# ScPnP 

This is an implementation of our published paper "ScPnP: A Non-iterative Scale Compensation Solution for PnP Problems".



## Usage:

In matlab: 

```
addpath ScPnPv;
addpath ScPnPz;
[ Rv,Tv ] = ScPnPv(vs,ps);
[ Rd,Td ] = ScPnPz(vs,ps);
```

where $v_s$ and $p_s$  are $3\times n$ matrices. $v_s$ are the measurement made by the camera and $p_s$ are the 3D points.



## Notes:

- This implementation is NOT efficient, reliable or versatile.  We are NOT responsible for any usage related to this code.
- In ScPnPv, there is a simple implementation of the symbolic computation.
- Our method obtains explicit expressions of the action matrix, which are quite long. Matlab takes a very long time to load these expressions. Fortunately, it seems that Matlab only loads them once.



## Improvement:

Comparing to our published paper, two improvements are made in this implementation:

### 1 

In our published paper, every unknown is solved from Dixion resultant respectively. However, the rest unknowns are already given in the eigenvectors as one computes the first unknown. Using eigenvector may save a bit computational time.

### 2

In our published paper, the high order variable, such as $q_0$，is set to 1 in the objective function. Lately we found that this operation can be postponed to significantly simplify the problem. Let $q_0$ still be a variable. The optimal condition reads:
$$
\frac{J}{\partial q_0}=0\\
\frac{J}{\partial q_1}=0\\
\frac{J}{\partial q_2}=0\\
\frac{J}{\partial q_3}=0\\
\frac{J}{\partial x_4}=0\\
$$

Due to scale ambiguity, the objective function does not change by scaling up $q_0,q_1,q_2,q_3$ with $\alpha$ and scaling down $x_4$ with $\frac{1}{\alpha^2}$. One can now substitute $q_0$ or $q_1,q_2,q_3$ with 1 and abandon one of the equation. The order of the top 4 equations is $3^{th}$ respectively while the order of the last one is $4^{th}$. So we abandon the last one.  

remark: Substituting $q_0$ in the objective function is equivalent to abandon the first equation.

As a result, the size of the Dixon matrix is reduced from $22\times22$ to $17\times17$. The size of the action matrix is reduced from $131\times 131$ to $91\times91$. The expressions are also much simpler. 

Most importantly, it requires roughly 60% computational time compare to our original method.

 

## email: 

xww@cad.zju.edu.cn

or

mczhj-whdx@126.com
