

<br/>
<br/>
<br/>

$$
b_{i}^{k}=
 \begin{cases}
   w^{k}  &\mbox{$i=o(k)$}\\
   -w^{k}  &\mbox{$i=d(k)$}\\
   0&\mbox{otherwise}
   \end{cases}
$$


***Minimize***
$\sum_{\tau\in \widetilde{\theta}} \sum_{(i,j)\in A} v_{ij}r_{ij}^{\tau}z_{\tau}+ \sum_{p\in P}\sum_{l\in L}\sum_{\tau\in \widetilde{\theta} _{pl}}  f_{pl}z _{\tau}+\sum_{k\in K} \sum_{(i,j)\in A} c_{ij}x_{ij}^{k}$

***Subject to:***
   $\sum_{j:(i,j)\in A} x_{ij}^{k}-\sum_{j:(j,i)\in A} x_{ji}^{k}=b_{i}^{k}
    \qquad \forall k\in K, \forall i\in N\qquad (1)$
   
$\sum_{k\in K}x_{ij}^{k} \leqslant \sum_{p\in P}\sum_{\tau\in \widetilde{\theta} _{p}} C_{p}r_{ij}^{\tau}z_{\tau} \qquad \forall(i,j)\in E\qquad(2)$

   $\sum_{\tau\in \widetilde{\theta}_{pl}}z_{\tau} \leqslant ub_{pl}\qquad \forall p\in P,\forall l\in L \qquad(3)$

   $x_{ij}^{k}\geqslant 0\qquad \forall k\in K,(i,j) \in A\qquad(4)$

$z_{\tau}\geqslant 0\qquad \forall \tau\in \widetilde{\theta}\qquad(5b)$

$\sum_{\tau\in \widetilde{\theta}}g_{s_{j}}^{\tau}z_{\tau} \leqslant U^{j} \qquad \forall j\in G^{u-}_{1}\qquad(6a)$

$\sum_{\tau\in \widetilde{\theta}}g_{s_{j}}^{\tau}z_{\tau} \geqslant L^{j} \qquad \forall j\in G^{u+}_{1}\qquad(6b)$


$\sum_{\tau\in \widetilde{\theta}_{p_{j} l_{j}}}g_{s_{j}}^{\tau}z_{\tau} \leqslant U^{j} \qquad \forall j\in G^{u-}_{2}\qquad(7a)$

$\sum_{\tau\in \widetilde{\theta}_{p_{j} l_{j}}}g_{s_{j}}^{\tau}z_{\tau} \geqslant L^{j} \qquad \forall j\in G^{u+}_{2}\qquad(7b)$

$\sum_{\tau\in \widetilde{\theta}_{p_{j} l_{j}}}r_{a_{j}}^{\tau}z_{\tau} \leqslant U^{j} \qquad \forall j\in G^{u-}_{3}\qquad(8a)$

$\sum_{\tau\in \widetilde{\theta}_{p_{j} l_{j}}}r_{a_{j}}^{\tau}z_{\tau} \geqslant L^{j} \qquad \forall j\in G^{u+}_{3}\qquad(8b)$

$x_{a_{j}}^{k_{j}} \leqslant w^{k_{j}}\sum_{\tau\in \widetilde{\theta}}r_{a_{j}}^{\tau}z_{\tau} \qquad \forall j\in cutset\qquad(9)$