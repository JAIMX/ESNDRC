Model
===================





Parameters
-------------
**卡车相关参数**

 - $\alpha:$ fixed cost per day per truck
 - $\beta$ :   transportation cost per package per unit distance
 - $s$ :  a sequence of index of different capacities
 - $C_{s}$:   capacity of $s th$ type of truck
 - $L:$ max number of legs allowed to be traveled by a truck
 -  $D:$ max distance allowed to be traveled by a truck
 -  $Speed:$ average speed of trucks, if necessary it can be truck specific
 -  $DrivingTimePerDay:$ driving time per day allowed for trucks
 - $b_{ij}^{\tau}=1$ if  $\tau$ contains arc(i,j)
 - $L_{so}:$  the number of trucks available starting from origin o with capacity of  $C_{s}$


**节点相关参数**

 - $q^{p}$: quantity of pickup and delivery demand $p $
 - $l_{i,j}$: distance of arc$(i,j) $

***Auxiliary graph***  $\qquad G^{'}(V^{'},A^{'})$

 - $V^{'}:$ for each $u \in V $, associate T vertices: $u_{1},u_{2},...,u_{T}$
 - $A_{T}=\left \{ (u_{t},u_{(t+1)mod T})|u\in  V , t\in\left \{ 1,2,...,T \right \}\right \}$
 - $\tilde{A}=\left \{ (u_{t},(w_{(t+t(u,w))mod T}))|(u,w)\in  A , t\in\left \{ 1,2,...,T \right \}\right \}$
 - $A^{'}=A_{T}\cup \tilde{A}$
 - cost: $A_{T}=0  \qquad\tilde{A}=l(u,w)$     
 

<br/>
<br/>
<br/>
<br/>

 
Decision variables
-------------
- $z_{\tau}$: the number of trucks choose to run in cycle $\tau$
- $x_{i,j}^{p}$: a split of demand $q^{p}$ shipped on arc $(i,j)\in \tilde{A}\cup A_{T}$


Sets
-------------
- $V$:  set of nodes
- $A$:  set of arcs
- $P$:  set of demand O-D pairs
- $S$:  set of index of different capacities of trucks

Indices
-------------
- $i,j$:  index of nodes
- $(i,j)$:  index of arcs
- $p$:  index of O-D pairs

Const
-------------


$$
b_{i}^{p}=
 \begin{cases}
   q^{p}  &\mbox{$i=o(p)$}\\
   -q^{p}  &\mbox{$i=d(p)$}\\
   0&\mbox{otherwise}
   \end{cases}
$$

<br/>
<br/>
<br/>
<br/>

***Minimize***
$\sum_{\tau\in \theta} \sum_{(i,j)\in \tilde{A}\cup A_{T}} \frac{\alpha l_{ij}b_{ij}^{\tau}z_{\tau}}{Speed*DrivingTimePerDay}+ \sum_{s\in S}\sum_{\tau\in \theta _{s}} \gamma _{s}z _{\tau}+\sum_{p\in P} \sum_{(i,j)\in \tilde{A}\cup A_{T}}\beta l_{ij}x_{ij}^{p}$

***Subject to:***
  
   
   $\sum_{(i,j)\in \delta ^{+}(i)} x_{ij}^{p}-\sum_{(j,i)\in \delta ^{-}(i)} x_{ji}^{p}=b_{i}^{p} \qquad \forall p\in P, i\in V^{'}\qquad (1)$
   
$\sum_{p\in P}x_{ij}^{p} \leqslant \sum_{s\in S}\sum_{\tau\in \theta _{s}} C_{s}b_{ij}^{\tau}z_{\tau} \qquad \forall(i,j)\in \tilde{A}\qquad(2)$

$\sum_{\tau\in \theta}b_{ij}^{\tau}z_{\tau} \leqslant1 \qquad\forall(i,j)\in \tilde{A}\qquad(3)$

   $\sum_{\tau\in \theta_{so}}z_{\tau} \leqslant L_{so}\qquad \forall s\in S,\forall o\in O \qquad(4)$
   

   
   $x_{ij}^{p}\geqslant 0\qquad \forall p\in P,(i,j) \in \tilde{A}\cup A_{T}\qquad(5)$
   $z_{\tau}\in Z\qquad \forall \tau\in \theta\qquad(6)$
   