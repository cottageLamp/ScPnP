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
- In our published paper, every unknown is solved from their Dixion resultant respectively.  Lately we found that, as computing the first unknown,  the rest of the unknowns are already obtained in the eigenvectors. The two methods return the same result but using eigenvector may save a bit computational time. In this code, we use eigenvectors in ScPnPv while keep ScPnPz compliance to our paper.
- In ScPnPv, there is a simple implementation of the symbolic computation.
- Our method obtains explicit expressions of the action matrix. Those expressions are quite long and take a long time for Matlab to load. Fortunately, it seems that Matlab only needs to load them once.



## email: 

xww@cad.zju.edu.cn

or

mczhj-whdx@126.com