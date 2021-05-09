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

where $v_s$ and $p_s$  are $3\times n$ matrices. $v_s$ contains the measurement seen by the camera and $p_s$ has the coordinates of the reference points.



## Notes:

- This implementation is NOT efficient, reliable or versatile.  We are NOT responsible for any usage related to this code.
- In ScPnPv, there is a simple implementation of the symbolic computation.
- Our method obtains explicit expressions of the action matrix. Those expressions are quite long and take a long time for Matlab to load. Fortunately, it seems that Matlab only loads them once.



## Improvement:

Comparing to our published paper, this implementation has made two improvements:

### 1 

In our published paper, every unknown is solved from their Dixion resultant respectively.  Lately we found that, as computing the first unknown,  the rest of the unknowns are already obtained in the eigenvectors. The two methods return the same result. However using eigenvector may save a bit computational time. In this code, we use eigenvectors in ScPnPv while keep ScPnPz compliance to our paper.

### 2

In our paper, the high order variable, such as $q_0$，is set to 1 in the objective function. Lately we found that this operation can be postponed and significantly simplify the problem. Let $q_0$ still be a variable. The optimal conditions reads:
$$
\frac{J}{\partial q_0}=0\\
\frac{J}{\partial q_1}=0\\
\frac{J}{\partial q_2}=0\\
\frac{J}{\partial q_3}=0\\
\frac{J}{\partial x_4}=0\\
$$
Due to scale ambiguity, the objective function keeps unchanged while one scales up $q_0,q_1,q_2,q_3$ with $\alpha$ and scales down $x_4$ with $\frac{1}{\alpha^2}$. One can now substitute $q_0$ or $q_1,q_2,q_3$ with 1 and abandon one of the equation. The order of the top 4 equations is $3^{th}$ while the order of the last one is $4^{th}$. So we abandon the last one.  

remark: Substituting $q_0$ in the objective function is equivalent to abandon the first equation.

As a result, the size of the Dixon matrix drops from $22\times22$ to $17\times17$. The size of the action matrix drops from $131\times 131$ to $91\times91$. The expressions are also much simpler. 

Most importantly, it only requires roughly 60% computational time compare to our original method to solve the PnP problem.

 

## email: 

xww@cad.zju.edu.cn

or

mczhj-whdx@126.com