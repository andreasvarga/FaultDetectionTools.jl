
# Synthesis paradigms 


The implemented computational procedures for the synthesis of fault detection filters share several computational paradigms, which are instrumental in developing generally applicable, numerically reliable and computationally efficient synthesis methods.
In what follows we shortly review these paradigms and discuss their roles in the synthesis procedures.

## Nullspace-based synthesis

An important synthesis paradigm is the use of the nullspace method as a first synthesis step to ensure the fulfillment of the decoupling conditions  $R_u(\lambda) = 0$ and $R_d(\lambda) = 0$. This can be done by choosing $Q(\lambda)$ of the form
```math
 Q(\lambda) = \overline Q_1(\lambda) Q_1(\lambda) , 
``` 
where the factor $Q_1(\lambda)$ is a left nullspace basis of the rational matrix
```math
G(\lambda) := \left[ \begin{array}{cc} G_u(\lambda) & G_d(\lambda) \\ I_{m_u} & 0 \end{array}\right] \, .
```
It follows
```math
[\,R_u(\lambda) \; R_d(\lambda)\,] = Q(\lambda)G(\lambda)= 0 \, .
```

The residual generator filter can be rewritten in the alternative form
```math
{\mathbf{r}}(\lambda) = \overline Q_1(\lambda)Q_1(\lambda)\left[ \begin{array}{c}
{\mathbf{y}}(\lambda)\\{\mathbf{u}}(\lambda)\end{array}\right] = \overline Q_1(\lambda) \overline{\mathbf{y}}(\lambda) \;, 
```
where
```math
\overline{\mathbf{y}}(\lambda) := Q_1(\lambda)\left[\begin{array}{c}
{\mathbf{y}}(\lambda)\\{\mathbf{u}}(\lambda)\end{array}\right] = \overline G_f(\lambda){\mathbf{f}}(\lambda) + \overline G_w(\lambda){\mathbf{w}}(\lambda)  \,,
```
with
```math
[\, \overline G_f(\lambda) \; \overline G_w(\lambda) \,] := Q_1(\lambda)
\left[ \begin{array}{cc} G_f(\lambda) & G_w(\lambda) \\ 0 & 0 \end{array}\right]\, .
```
With this first preprocessing step, the original problems formulated for a system with control and disturbance inputs can be reformulated for the above reduced system (without control and disturbance inputs),  for which we have to determine the TFM $\overline Q_1(\lambda)$ of the simpler fault detection filter. For the details of the implemented computational approach see Section 7.4 of [1]. 

For the computation of nullspace bases, functions available in the [DescriptorSystems](https://github.com/andreasvarga/DescriptorSystems.jl) package are employed. 



## Using filter updating techniques

In all implemented synthesis procedures, the TFM of the resulting filter $Q(\lambda)$ can be expressed in a factored form as
```math
Q(\lambda) = Q_K(\lambda) \cdots Q_2(\lambda)Q_1(\lambda) \, , 
```
where $Q_1(\lambda)$ is a left nullspace basis of the above defined $G(\lambda)$, satisfying $Q_1(\lambda)G(\lambda) = 0$, and  $Q_1(\lambda)$, $Q_2(\lambda)Q_1(\lambda)$, $\ldots$, can be interpreted as partial syntheses addressing specific requirements. Since each partial synthesis may represent a valid fault detection filter, this approach can be flexibly used  for employing or combining different synthesis techniques.

The determination of $Q(\lambda)$ in the above factored form can be formulated as a $K$-step synthesis procedure based on successive updating of an initial filter $Q(\lambda) = Q_1(\lambda)$ and the nonzero terms of its corresponding  internal form
```math
R(\lambda) := [\, R_f(\lambda) \; R_w(\lambda)\,] = Q_1(\lambda)\left[ \begin{array}{cc} G_f(\lambda) & G_w(\lambda) \\  0 & 0 \end{array}\right]
```
as follows: for $k = 2, \ldots, K$, determine $Q_k(\lambda)$ using the current $Q(\lambda)$ and $R(\lambda)$ and  then perform the updating as $Q(\lambda) \leftarrow Q_k(\lambda)Q(\lambda)$ and 
$R(\lambda) \leftarrow Q_k(\lambda)R(\lambda)$. These updating operations are efficiently performed using state-space description based formulas.
The main benefit of using explicit state-space based updating formulas is the possibility to ensure at each step the cancellation of a maximum number of poles and zeros between the factors. This allows to keep the final order of the filter $Q(\lambda)$ as low as possible. See Section 7.3 of [1] for a discussion of additional aspects. 

## Least order synthesis based on admissibility conditions

The least order synthesis of fault detection filters means to determine filters $Q(\lambda)$ with the least possible orders, to help in reducing the computational burden associated with their real-time implementations. The main tool to achieve least order synthesis is the solution of suitable _minimal cover problems_ (see Section 7.5 of [1]). If $X_1(\lambda)$ and $X_2(\lambda)$ are rational matrices of the same column dimension,  then the _left minimal cover problem_ is to find $X(\lambda)$ and $Y(\lambda)$ such that
```math
X(\lambda) = X_1(\lambda) + Y(\lambda) X_2(\lambda) , 
```
and the order of $[\,X(\lambda) \; Y(\lambda) \,]$ is minimal.


A typical second step in many synthesis procedures is to choose $Q_2(\lambda)$ such that
the product $Q_2(\lambda)Q(\lambda)$ has least dynamical order and, simultaneously, a certain _admissibility_ condition is fulfilled (usually involving the nonzero TFMs $R_f(\lambda)$ and $R_w(\lambda)$). For example, for the solution of the AFDP, the rows of $Q(\lambda) := Q_1(\lambda)$ form a left nullspace basis and the employed admissibility conditions  are
```math
Q_2(\lambda)R_{f_j}(\lambda) \not = 0, \;\;j = 1, \ldots, m_f, 
```
which  guarantee the detectability of all fault components, and additionally, the full row rank requirement on $Q_2(\lambda)R_w(\lambda)$. The latter condition is only imposed for convenience, to simplify the subsequent computational steps. 

The determination of candidate solutions $Q_2(\lambda)$ such that $Q_2(\lambda)Q(\lambda)$ has least order  can be done by solving left minimal cover problems, where $X_1(\lambda)$ and $X_2(\lambda)$ represent disjoint subsets of basis vectors, such that: $Q(\lambda) = \left[\begin{smallmatrix} X_1(\lambda) \\ X_2(\lambda) \end{smallmatrix}\right]$, $Q_2(\lambda) = [\, I \;\; Y(\lambda)\,]$,  and $X(\lambda) = Q_2(\lambda)Q(\lambda)$ and $Y(\lambda)$ represent the solution of the left cover problem. A systematic search over increasing orders of candidate solutions can be performed and the search stops when the admissibility conditions are fulfilled.

State-space representation based computational methods for the solution of minimum dynamic cover problems are described in Sections 7.5 and 10.4 of [1], together with explicit updating formulas of the state-space realizations of  $Q(\lambda)$ and $R(\lambda)$.

For the solution of minimal cover problems, functions available in the [DescriptorSystems](https://github.com/andreasvarga/DescriptorSystems.jl) package are employed. 




## Stabilization via left coprime factorization

A desired dynamics of the resulting final filters $Q(\lambda)$ and $R(\lambda)$ can be enforced by choosing a suitable invertible factor $M(\lambda)$, such that  $M(\lambda)[\, Q(\lambda) \; R(\lambda) \,]$ has desired poles. This can be achieved by computing a left coprime factorization
```math
[\, Q(\lambda) \; R(\lambda) \,] = M^{-1}(\lambda) [\, N_Q(\lambda) \; N_R(\lambda) \,] 
```
with $M(\lambda)$ and $[\, N_Q(\lambda) \; N_R(\lambda) \,]$ coprime and having arbitrary stable poles, and then performing the updating operations $Q(\lambda) \leftarrow N_Q(\lambda)$ and 
$R(\lambda) \leftarrow N_R(\lambda)$. The stabilization via a left coprime factorization is usually performed as the last step of the synthesis procedures. 

For a detailed discussion of the coprime factorization based stabilization approach see Sections 7.6 and 10.4 of [1], where coprime factorization methods, based on recursive pole dislocation techniques are described, which produce directly the numerator factors $N_Q(\lambda)$ and $N_R(\lambda)$ (thus implicitly perform  the updating operations as well).


For the solution of coprime factorizations, functions available in the [DescriptorSystems](https://github.com/andreasvarga/DescriptorSystems.jl) package are employed. 



### References

[1]   A. Varga, Solving Fault Diagnosis Problems – Linear Synthesis Techniques, Vol. 84 of Studies in Systems, Decision and Control, Springer International Publishing, 2017.




