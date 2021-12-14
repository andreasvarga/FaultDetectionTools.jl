
# Fault detection and diagnosis basics 

## Synthesis models

The plant models underlying all implemented synthesis methods (also called _synthesis models_) are linear time-invariant system models, where the faults
are equated with special (unknown) disturbance inputs. 
An important class of models with additive
faults arises when defining the fault signals for two
main categories of faults, namely, actuator and sensor
faults. Two basic forms of synthesis models are used.

The input-output plant model with additive faults  has the form
```math
{\mathbf{y}}(\lambda) = G_u(\lambda){\mathbf{u}}(\lambda) +
G_d(\lambda){\mathbf{d}}(\lambda) + G_f(\lambda){\mathbf{f}}(\lambda) +
G_w(\lambda){\mathbf{w}}(\lambda) ,
```
where  ``{\mathbf{y}}(\lambda)``, ``{\mathbf{u}}(\lambda)``,
``{\mathbf{d}}(\lambda)``, ``{\mathbf{f}}(\lambda)``, and
``{\mathbf{w}}(\lambda)`` (with boldface notation), are the Laplace-transformed (in the continuous-time case) or Z-transformed (in the discrete-time case)  ``p``-dimensional system output vector ``y(t)``,
``m_u``-dimensional control input vector ``u(t)``,
``m_d``-dimensional disturbance vector ``d(t)``,
``m_f``-dimensional fault vector ``f(t)`` and
``m_w``-dimensional noise vector ``w(t)``,
respectively, and where ``G_u(\lambda)``, ``G_d(\lambda)``,
``G_f(\lambda)``   and ``G_w(\lambda)`` are proper _transfer function matrices_ (TFMs) from the respective
inputs to outputs.  Input-output models with additive faults of the above form are useful in formulating various fault diagnosis problems, in deriving general solvability conditions and in describing conceptual synthesis procedures. However, these models are generally not suited for numerical computations, due to the  potentially high sensitivity of polynomial-based model representations.

For computational purposes, instead of the above input-output model with the compound TFM $[\, G_u(\lambda) \; G_d(\lambda) \; G_f(\lambda) \; G_w(\lambda) \,]$, an equivalent state-space model is used having the form
```math
\begin{array}{rcl} E\lambda x(t) &=& Ax(t) + B_u u(t) + B_d d(t) + B_f f(t) + B_w w(t) \, ,  \\
y(t) &=&Cx(t) + D_u u(t) + D_d d(t) +  D_f f(t) + D_w w(t)  \, ,
\end{array}
```
 with the $n$-dimensional state vector $x(t)$, where $\lambda
x(t) := \dot{x}(t)$ or $\lambda x(t) := x(t+1)$ depending on the
type of the system, continuous- or discrete-time, respectively.  The matrix $E$ is generally invertible and is frequently  taken as $E = I_n$. Plant models of the above state-space form often arise from the linearization of nonlinear dynamic plant models in specific operation points and for  fixed values of plant parameters. The noise inputs generally  account  for the effects of uncertainties (e.g., inherent  variabilities in operating points and parameters).

To indicate the input-output equivalence of the two types of models, we use the notation
```math
[\, G_u(\lambda) \; G_d(\lambda) \; G_f(\lambda) \; G_w(\lambda) \,] = \left[ \begin{array}{c|cccc} A-\lambda E & B_u & B_d & B_f & B_w \\  \hline
C & D_u & D_d & D_f & D_w \end{array}\right] ,
```
which stays for the following relations between the elements of the two representations:
```math
\begin{array}{lll} G_u(\lambda) &=& C(\lambda E-A)^{-1}B_u+D_u \\
G_d(\lambda) &=& C(\lambda E-A)^{-1}B_d+D_d \\
G_f(\lambda) &=& C(\lambda E-A)^{-1}B_f + D_f\\
G_w(\lambda) &=& C(\lambda E-A)^{-1}B_w + D_w 
\end{array}
```

The state-space synthesis model underlies the definition of the [`FDIModel`](@ref) object used in the
functions for the synthesis of fault diagnosis filters.

## Residual generation

A nonzero fault signal $f \neq 0$ signifies a deviation from the normal behaviour of the plant
due to an unexpected event (e.g., physical component failure or
supply breakdown). Generally, the occurrence of a fault
must be detected as early as possible to prevent further degradation
of the plant behaviour. _Fault detection and diagnosis_ (FDD) is concerned with one or more of the following aspects: the detection of the occurrence of any fault (_fault detection_), the
localization of detected faults (_fault isolation_), the reconstruction of the fault signal (_fault estimation_) and the classification of the
detected faults and determination of their characteristics (_fault identification_). The later aspect is not addressed in this package. 

A FDD system is a device (usually based on a collection of real-time processing algorithms)
suitably set-up to fulfill the above tasks. 
The main component of any FDD system is the _residual
generator_ (or _fault detection filter_),
which produces residual signals grouped in a $q$-dimensional vector $r$ by processing the available
measurements $y$ and the known values of control inputs $u$. The role of the residual signals
is to indicate the presence  or absence of faults, and therefore the residual $r$ must be equal (or
close) to zero in the absence of faults and significantly different
from zero after a fault occurs.

For decision-making when solving fault detection problems, a suitable measure of the residual magnitude is generated in a scalar _evaluation signal_ $\theta$  (e.g., $\theta = \|r\|$), which is then used to set a _decision variable_, say  $\iota$, as follows: $\iota= 1$, if $\theta > \tau$ for a detected fault, and $\iota= 0$ if $\theta \leq \tau$,  for the lack of faults, where $\tau$ is a given detection threshold.

For decision-making when solving faul isolation problems, $r(t)$ is generally a structured vector with, say $n_b$ components $r^{(i)}(t)$, $i = 1, \ldots, n_b$, and  $\theta$ and $\iota$ are $n_b$-dimensional vectors, with $\theta_i$ representing a measure of the magnitude of the $i$-th residual component (e.g., $\theta_i = -\|r^{(i)}\|$). The $i$-th component of the binary signature vector $\iota$ is set $\iota_i = 1$ or $\iota_i = 0$ corresponding to a fired (i.e, $\theta_i > \tau$) or not fired (i.e., $\theta_i \leq \tau$)  component $r^{(i)}(t)$, respectively.

A linear residual generator employed in a FDD system  has the input-output form
```math
{\mathbf{r}}(\lambda) = Q(\lambda)\left[ \begin{array}{c}
{\mathbf{y}}(\lambda)\\{\mathbf{u}}(\lambda)\end{array}\right] = 
Q_y(\lambda){\mathbf{y}}(\lambda) + Q_u(\lambda){\mathbf{u}}(\lambda) ,
```
where
$Q(\lambda) = [Q_y(\lambda) \; Q_u(\lambda)]$ is the TFM of the filter. For a physically
realizable filter,  $Q(\lambda)$ must be _stable_ (i.e., with all its
poles having negative real parts for a continuous-time system
or magnitudes less than one for a discrete-time system). The
_order_ of $Q(\lambda)$ is the dimension of the state vector
of a minimal state-space realization of $Q(\lambda)$. The
dimension $q$ of the residual vector $r(t)$ depends on the
fault diagnosis problem to be solved. The above input-output representation is called the _implementation form_ of the fault detection filter and is the basis of its real-time implementation.



The residual signal $r(t)$ generally depends
via the system outputs $y(t)$ of all system inputs $u(t)$,
$d(t)$, $f(t)$ and $w(t)$. The _internal form_ of the filter is obtained by replacing in the above equation ${\mathbf{y}}(\lambda)$ by its expression in the synthesis model, and is given by
```math
{\mathbf{r}}(\lambda) = R(\lambda)\left[ \begin{array}{c}{\mathbf{u}}(\lambda)\\
{\mathbf{d}}(\lambda) \\ {\mathbf{f}}(\lambda) \\ {\mathbf{w}}(\lambda)\end{array}\right] =  
R_u(\lambda){\mathbf{u}}(\lambda) +
R_d(\lambda){\mathbf{d}}(\lambda) + 
R_f(\lambda){\mathbf{f}}(\lambda) + R_w(\lambda){\mathbf{w}}(\lambda)
,
```
where
```math
R(\lambda) = [\, R_u(\lambda) \mid  R_d(\lambda) \mid  R_f(\lambda) \mid  R_w(\lambda)\,] := 
Q(\lambda)  \left[ \begin{array}{c|c|c|c} G_u(\lambda) & G_d(\lambda) & G_f(\lambda)  & G_w(\lambda) \\
         I_{m_u} & 0 & 0 & 0 \end{array}\right] .
```
For a successfully designed filter $Q(\lambda)$, all TFMs in the
corresponding internal form $R(\lambda)$ are  stable, and  fulfil  specific fault diagnosis requirements.


The basic functionality  of a well-designed fault detection filter is to ensure the lack of false alarms, in the case when no faults occurred, and the lack of missed detection of faults, in the case of occurrence of a fault. The first requirement is fulfilled if, in the presence of noise, the signal norm $\|r\|$ is  sufficiently small for all possible control, disturbance and noise inputs. The requirement on the lack of missed detections is fulfilled provided $\|r\|$ is sufficiently large for any fault of sufficiently large magnitude for all possible control, disturbance and noise inputs.

The fault detection filter in the implementation form underlies the definition of the [`FDFilter`](@ref) object, while the internal form underlies the definition of the [`FDFilterIF`](@ref) object. These objects are generated by several synthesis functions. 

For the isolation of faults, a bank of residual generator filters is employed which is formed by stacking a bank of $n_b$ filters of the form
```math
{\mathbf{r}}^{(i)}(\lambda) = Q^{(i)}(\lambda)\left[ \begin{array}{c}
{\mathbf{y}}(\lambda)\\{\mathbf{u}}(\lambda)\end{array}\right] \, ,
```
where the $i$-th filter $Q^{(i)}(\lambda)$ generates the corresponding $i$-th residual component 
$r^{(i)}(t)$ (scalar or vector). The internal form of the $i$-th filter $Q^{(i)}(\lambda)$ is 
```math
{\mathbf{r}}^{(i)}(\lambda) = R^{(i)}(\lambda)\left[ \begin{array}{c}{\mathbf{u}}(\lambda)\\
{\mathbf{d}}(\lambda) \\ {\mathbf{f}}(\lambda) \\ {\mathbf{w}}(\lambda)\end{array}\right] =  
R^{(i)}_u(\lambda){\mathbf{u}}(\lambda) +
R^{(i)}_d(\lambda){\mathbf{d}}(\lambda) + 
R^{(i)}_f(\lambda){\mathbf{f}}(\lambda) + R^{(i)}_w(\lambda){\mathbf{w}}(\lambda)
,
```
where
```math
R^{(i)}(\lambda) = [\, R^{(i)}_u(\lambda) \mid  R^{(i)}_d(\lambda) \mid  R^{(i)}_f(\lambda) \mid  R^{(i)}_w(\lambda)\,] := 
Q^{(i)}(\lambda)  \left[ \begin{array}{c|c|c|c} G_u(\lambda) & G_d(\lambda) & G_f(\lambda)  & G_w(\lambda) \\
         I_{m_u} & 0 & 0 & 0 \end{array}\right] .
```

This leads to the following structured residual vector $r(t)$ and block-structured filters $Q(\lambda)$ and $R(\lambda)$
```math
r(t) = \left[ \begin{array}{c} r^{(1)}(t)\\ \vdots \\ r^{(n_b)}(t) \end{array}\right] , \;
Q(\lambda) = \left[ \begin{array}{c} Q^{(1)}(\lambda)\\ \vdots \\ Q^{(n_b)}(\lambda) \end{array}\right]  , \; R(\lambda) = \left[ \begin{array}{c} R^{(1)}(\lambda)\\ \vdots \\ R^{(n_b)}(\lambda) \end{array} \right]  \,.
```

The above bank of fault detection filters is the basis of the definition of the fault detection and isolation object [`FDIFilter`](@ref) and its internal form [`FDIFilterIF`](@ref).

# Structure matrix

Consider $R_f(\lambda)$, the TFM from the fault inputs $f$ to residual $r$ in the internal form, and assume $R_f(\lambda)$ is an $n_b\times m_f$ block-structured TFM of the form
```math
R_f(\lambda) = \left[ \begin{array}{ccc} R^{(1)}_{f_1}(\lambda)& \cdots &R^{(1)}_{f_{m_f}}(\lambda) \\
\vdots & \ddots & \vdots \\
 R^{(n_b)}_{f_1}(\lambda)& \cdots &R^{(n_b)}_{f_{m_f}}(\lambda) \end{array}\right] \, , 
```
where the $(i,j)$-th block of $R_f(\lambda)$ is defined as 
```math
R^{(i)}_{f_j}(\lambda) := Q^{(i)}(\lambda) \left[ \begin{array}{c} G_{f_j}(\lambda) \\ 0 \end{array}\right] .
```
Here, $Q^{(i)}(\lambda)$ is either the $i$-th row of the filter $Q(\lambda)$, in which case $R^{(i)}_{f_j}(\lambda)$ is the $(i,j)$-th element of $R_f(\lambda)$, or the $i$-th block row of $Q(\lambda)$ corresponding to the $i$-th filter in a bank of $n_b$ filters. 
In both cases, $R^{(i)}_{f_j}(\lambda)$ describes how the $j$-th fault $f_j$ influences the $i$-th (scalar or vector) component of the residual $r$.

We associate to the block structured $R_f(\lambda)$ the $n_b\times m_f$ binary
_structure matrix_ $S_{R_f}$, whose $(i,j)$-th element is defined as
```math
\begin{array}{llrll} S_{R_f}(i,j) &=& 1 & \text{ if } & R^{(i)}_{f_j}(\lambda) \not=0 \; ,\\
S_{R_f}(i,j) &=& 0 & \text{ if } & R^{(i)}_{f_j}(\lambda) =0 \, .
\end{array}
```
If $S_{R_f}(i,j) = 1$
then we say that the $i$-th residual component is sensitive to the $j$-th fault $f_j$, while if $S_{R_f}(i,j) = 0$ then the $j$-th fault $f_j$  is decoupled from $i$-th residual component. The $m_f$ columns of $S_{R_f}$ are called _fault signatures_ and play a crucial role in the decision making 
for fault isolation. Since each nonzero column of $S_{R_f}$ is associated with the corresponding fault input, fault isolation can be performed by comparing the resulting binary decision vector $\iota$ in the FDD system (i.e., the signatures of fired or not fired residual components) with the fault signatures coded in the columns of $S_{R_f}$. The rows of $S_{R_f}$ play an important role in solving FDI synthesis problems and are called _specifications_. 

The above definition of the structure matrix $S_{R_f}$ is associated with the zero/nonzero blocks of the TFM $R_f(\lambda)$ and is also known as the _weak structure matrix_. 
The _strong structure matrix_ is related to the zero/nonzero blocks of the frequency response of $R_f(\lambda)$ evaluated for a set of 
relevant complex frequencies $\Omega$ characterizing the classes of persistent fault inputs. For example, to a _real_ frequency $\omega$ which characterizes sinusoidal faults, the corresponding _complex_ frequency in $\Omega$ is $j\omega$ for a continuous-time system or $\exp(j\omega T_s)$ for a discrete-time system with sampling time $T_s$ (thus, the null frequency characterizes constant faults). The strong structure matrix is defined as
```math
\begin{array}{llrll} S_{R_f}(i,j) &=& 1 & \text{ if } & R^{(i)}_{f_j}(\lambda_z) \not=0 \;  \text{ for all }  \lambda_z \in \Omega ,\\
S_{R_f}(i,j) &=& 0 & \text{ if } & R^{(i)}_{f_j}(\lambda_z) = 0 \, \text{ for any } \lambda_z \in \Omega .
\end{array} 
```

For the determination of the weak or strong structure matrix 
the function [`fditspec`](@ref) is available. Alternatively, the function [`fdisspec`](@ref) can be used to determine the strong structure matrix. 

When solving fault isolation problems, the choice of the desired structure matrix $S_{R_f}$ is an important aspect. The function [`fdigenspec`](@ref) allows to compute the maximally achievable structure matrix for a given synthesis model. The function [`fdichkspec`](@ref) can be employed to check the feasibility of a set of FDI specifications.


# Fault diagnosis problems

Six basic fault diagnosis problems are formulated in what follows and their solutions are addressed
by the implemented synthesis functions. To fulfill the basic requirement for the lack of false alarms
in the presence of arbitrary control and disturbance inputs, for all problems we require that by a suitable choice of a stable fault detection filter $Q(\lambda)$, we achieve that  the residual signal $r(t)$ is fully decoupled from the control input $u(t)$ and disturbance input $d(t)$. Thus, the following _decoupling conditions_ must be generally fulfilled:
```math
\begin{array}{ll}
  (i) & R_u(\lambda) = 0 ,\\
  (ii) & R_d(\lambda) = 0 .
\end{array}
```
Since the effect  of a nonzero  noise input $w(t)$ can usually not be fully decoupled from the residual $r(t)$, an additional requirement is that the influence of the noise signal $w(t)$ is negligible. Thus, the following _noise attenuation condition_ has to be also fulfilled:
```math
\begin{array}{ll}
 (iv) & R_w(\lambda) \approx 0, \;\; \textrm{with} \; R_w(\lambda) \; \textrm{stable.}
\end{array} 
```

The condition $R_w(\lambda) \approx 0$ expresses the requirement  that the transfer gain $\|R_w(\lambda)\|$ (measured by any suitable norm) can be made arbitrarily small and is intended to avoid missed detections in the presence of noise inputs. 

For particular fault diagnosis problems specific requirements on $R_f(\lambda)$ have to be additionally fulfilled.

We distinguish between exact and approximate solutions of fault diagnosis problems. The _exact_ problems impose no conditions regarding noise inputs, excepting the stability of $R_w(\lambda)$ in the case when $w \not\equiv 0$. The _approximate_ problems address the case of nonzero noise inputs
by employing special techniques to reduce their effects. For both cases we assume the following general internal form of the filter $Q(\lambda)$ 
```math 
r(\lambda) =
R_u(\lambda){\mathbf{u}}(\lambda) +
R_d(\lambda){\mathbf{d}}(\lambda)    + R_f(\lambda){\mathbf{f}}(\lambda) + R_w(\lambda){\mathbf{w}}(\lambda) 
```


### Exact fault detection problem (EFDP)
**EFDP:** Determine a stable residual generator $Q(\lambda)$ such that
```math
\begin{array}{rl} (i) & R_u(\lambda) = 0 \\ (ii) & R_d(\lambda) = 0 \\
  (iii) &R_{f_j}(\lambda) \not = 0,\; j = 1, \ldots, m_f \;\; \text{with} \;\; R_f(\lambda) \;\; {\color{magenta} \text{stable}} \\
  (iv) &R_w(\lambda) \;\; {\color{magenta} \text{stable}} 
  \end{array}
```
Condition $(iii)$ expresses the _complete fault detectability_ condition [1].

The EFDP can be formulated with the stronger requirement that the columns of $R_f(\lambda)$ do not vanish for a set of 
relevant complex frequencies $\Omega$ characterizing the classes of persistent fault inputs:

**EFDP with strong fault detectability:**  For a given set of complex frequencies $\Omega$, 
determine a stable	residual generator $Q(\lambda)$ such that
```math
\begin{array}{rl} (i) & R_u(\lambda) = 0 \\ (ii) & R_d(\lambda) = 0 \\
  (iii) &R_{f_j}(\lambda_z) \not = 0, \forall \lambda_z \in \Omega, \; j = 1, \ldots, m_f \;\; \text{with} \;\; R_f(\lambda) \;\; {\color{magenta} \text{stable}} \\
  (iv) &R_w(\lambda) \;\; {\color{magenta} \text{stable}} 
  \end{array}
```
Condition $(iii)$ expresses the _complete strong fault detectability_ condition [1].

When solving fault detection problems, it is important to assess the sensitivity of the residual signal to individual fault components. The assessment of  _complete fault detectability_ can be done by checking 
$\| R_{f}(\lambda) \|_{\infty -}  > 0$, where
```math
\| R_{f}(\lambda) \|_{\infty -} := \min_j \|R_{f_j}(\lambda)\|_\infty  
```
is the $\mathcal{H}_{\infty -}$-index defined in [1], as a measure of the degree of complete fault detectability. If $\| R_{f}(\lambda) \|_{\infty -}  = 0$, then at least one fault component is not detectable in the residual signal $r$. The assessment of the _strong complete fault detectability_ with respect to a set of frequencies contained in a set $\Omega$ comes down to check  $R_{f_j}(\lambda_z) \neq 0$, for $\forall \lambda_z \in \Omega$ and for $j = 1, \ldots, m_f$. Alternatively, the assessment of strong complete fault detectability can be done by checking 
$\| R_{f}(\lambda) \|_{\Omega -}  > 0$, where
```math
\| R_{f}(\lambda) \|_{\Omega -} := \min_{j}   \{ \inf_{\lambda_z \in \Omega} \|R_{f_j}(\lambda_z)\|_2 \}
```
is the (modified) $\mathcal{H}_{\infty -}$-index defined over the frequencies contained in $\Omega$ (see [1]). Since nonzero values of $\| R_{f}(\lambda) \|_{\infty -}$ or $\| R_{f}(\lambda) \|_{\Omega -}$ are not invariant to scaling (e.g., when replacing $Q(\lambda)$ by $\alpha Q(\lambda)$), these quantities are less appropriate to quantitatively assess the degrees of complete detectability.

The function [`fdhinfminus`](@ref) can be employed to evaluate $\|R_{f}(\lambda) \|_{\infty -} $ and $\| R_{f}(\lambda) \|_{\Omega -}$. 

A performance measure associated to a fault detection filter $Q(\lambda)$ which solves the EFDP is a scaling independent measure of the complete fault detectability called the _fault sensitivity condition_ and is defined (over all frequencies) as
```math
 J_1 =   \| R_{f}(\lambda) \|_{\infty -} / \max_j \|R_{f_j}(\lambda)\|_\infty. 
```
Similarly, scaling independent measure of the strong complete fault detectability 
can be defined over the frequencies contained in the set $\Omega$ as
```math
\widetilde J_1 =   \| R_{f}(\lambda) \|_{\Omega -} / \max_{j}\{ \sup_{\lambda_z \in \Omega} \|R_{f_j}(\lambda_z)\|_2 \}. 
```
We have that $0 < J_1 \leq 1$ and $0 < \widetilde J_1 \leq 1$ and a value of $J_1$ (or of $\widetilde J_1$) near to 1, indicates nearly equal sensitivities of residual to all fault components, and makes easier the choice of suitable thresholds for fault detection. On contrary, a small value of $J_1$ (or of $\widetilde J_1$) indicates potential difficulties in detecting some components of the fault vector, due to a very low sensitivity of the residual to these fault components. In such cases, employing fault detection filters with several outputs ($q > 1$) could be advantageous.

For the synthesis of fault detection filters which solve the EFDP the function [`efdsyn`](@ref) is available and for the evaluation of the  fault sensitivity condition the function [`fdiscond`](@ref) is available. 

### Exact fault detection and isolation problem (EFDIP)


**EFDIP:** Given a structure matrix $S$, determine a stable
	residual generator $Q(\lambda)$ such that
```math
\begin{array}{rl} (i) & R_u(\lambda) = 0 \\ (ii) & R_d(\lambda) = 0 \\
  (iii) & S_{R_{f}} = S \;\; \text{with} \;\; R_f(\lambda) \;\; {\color{magenta} \text{stable}} \\
  (iv) &R_w(\lambda) \;\; {\color{magenta} \text{stable}} 
  \end{array}
```
The condition (iii) expresses the $S$ _fault isolability
property_ [1].

The solution of the EFDIP can be addressed by solving $n_b$ suitably formulated EFDPs. The $i$-th EFDP arises by reformulating the $i$-th EFDIP for determining the $i$-th  filter $Q^{(i)}(\lambda)$ for a structure matrix which is the $i$-th row of $S$. This can be accomplished by redefining the fault components corresponding to zero entries in the $i$-th row of $S$ as additional disturbance inputs to be decoupled in the $i$-th residual component $r^{(i)}(t)$.

When solving fault detection and isolation problems with a targeted structure matrix $S$, we obtain partitioned filters (see above) and we can define for each individual filter an associated fault sensitivity condition number. Let $f^{(i)}$ be formed from the subset of faults corresponding to nonzero entries in the $i$-th row of $S$ and let $R_{f^{(i)}}^{(i)}(\lambda)$ be formed from the corresponding columns of $R_{f}^{(i)}(\lambda)$.  To characterize the  complete fault detectability of the subset of faults corresponding to nonzero entries in the $i$-th row of $S$ we can define the fault sensitivity condition number of the $i$-th filter as
```math
J_1^{(i)} =   \big\| R_{f^{(i)}}^{(i)}(\lambda) \big\|_{\infty -} / \max_j \big\|R_{f_j}^{(i)}(\lambda)\big\|_\infty . 
```
Similarly, to characterize the strong complete fault detectability of the subset of faults corresponding to nonzero entries in the $i$-th row of $S$,  we define the fault condition number of the $i$-th filter as
```math
\widetilde J_1^{(i)} =   \big\| R_{f^{(i)}}^{(i)}(\lambda) \big\|_{\Omega -} / \max_{j}\{ \sup_{\lambda_z \in \Omega} \big\|R_{f_j}^{(i)}(\lambda_z)\big\|_2 \} . 
```
For the synthesis of fault detection filters which solve the EFDIP the function [`efdisyn`](@ref) is available and for the evaluation of the  fault sensitivity condition the function [`fdiscond`](@ref) is available. 


### Approximate fault detection problem (AFDP)
**AFDP:** Determine a stable residual generator $Q(\lambda)$ such that
```math 
\begin{array}{rl} (i) & R_u(\lambda) = 0 \\ (ii) & R_d(\lambda) = 0 \\
  (iii) &R_{f_j}(\lambda) \not = 0,\; j = 1, \ldots, m_f \;\; \text{with} \;\; R_f(\lambda) \;\; {\color{magenta} \text{stable}} \\
  (iv) & \|R_w(\lambda) \| \approx 0 \;\; \text{with} \;\; R_w(\lambda) \;\; {\color{magenta} \text{stable}}
  \end{array}
```
The performance of a fault detection filter $Q(\lambda)$ which solves the AFDP can be characterized by the _fault-to-noise gap_ defined as 
```math 
J_2 := \|R_{f}(\lambda) \|_{\infty -} / \|R_w(\lambda)\|_{\infty}
```
By convention, $J_2 = 0$ if $\| R_{f}(\lambda) \|_{\infty -} = 0$ and $J_2 = \infty$ if $\| R_{f}(\lambda) \|_{\infty -} > 0$ and $\| R_{w}(\lambda) \|_{\infty} = 0$ (e.g., when solving exact synthesis problems without noise inputs). 
A finite frequency variant of the above criterion, which allows to address strong fault detectability aspects for a given set $\Omega$ of relevant frequencies is
```math  
\widetilde J_2 = \| R_{f}(\lambda) \|_{\Omega -} / \| R_{w}(\lambda) \|_{\infty} . 
```
The higher the value of $J_2$ (or $\widetilde J_2$), the easier is to choose suitable thresholds to be used for fault detection purposes in the presence of noise. Therefore, the maximization of the above gaps is a valuable goal in improving the fault detection capabilities of the fault diagnosis system in the presence of exogenous noise.


The value $\eta = J_2$ can be used to determine an estimation of the minimum size $\delta_{f,min}$ of detectable faults, as $\delta_{f,min} = \delta_w / \eta$, where, for $\delta_w$ is an upper bound
on the magnitude of the noise input $w(t)$ (i.e., $\|w \| \leq \delta_w$). 
The resulting value of $\delta_{f,min}$ can be used to assess the practical usefulness of any solution, and the maximization of the gap $J_2$ is always a meaningful goal for the synthesis of fault detection filters. Note that $J_2 = \infty$ for a filter $Q(\lambda)$ solving an EFDP with $w \equiv 0$. 

For the synthesis of fault detection filters which solve the AFDP the function [`afdsyn`](@ref) is available and for the evaluation of the  fault sensitivity condition the function [`fdif2ngap`](@ref) is available. 

### Approximate fault detection and isolation problem (AFDIP)

Let $S$ be a desired $n_b\times m_f$ structure matrix targeted to be achieved by using a structured fault detection filter $Q(\lambda)$ with $n_b$ row blocks (see also the formulation of the EFDIP) and let $R_f(\lambda)$ be the corresponding $n_b\times m_f$ block-structured TFM. $R_f(\lambda)$
can be additively decomposed as $R_f(\lambda) = \widetilde  R_f(\lambda) + \overline R_f(\lambda)$, where  $\widetilde  R_f(\lambda)$ and $\overline R_f(\lambda)$ have the same block structure as $R_f(\lambda)$ and have their $(i,j)$-th blocks defined as
```math 
\widetilde  R^{(i)}_{f_j}(\lambda) := S_{ij}R^{(i)}_{f_j}(\lambda), \quad \overline R^{(i)}_{f_j}(\lambda) := (1-S_{ij})R^{(i)}_{f_j}(\lambda) 
```
To address the approximate fault detection and isolation problem, we will target to enforce for
the part $\widetilde R_f(\lambda)$ of $R_f(\lambda)$ the desired structure matrix $S$, while the part $\overline R_f(\lambda)$ must (ideally) be negligible.

**AFDIP:** Given a structure matrix $S$, determine a stable residual generator $Q(\lambda)$ such that
```math 
\begin{array}{rl} (i) & R_u(\lambda) = 0 \\ (ii) & R_d(\lambda) = 0 \\
  (iii) &  S_{\widetilde R_f} = S, \; \overline R_f(\lambda) \approx 0\;\; \text{with} \;\; R_f(\lambda) \;\; {\color{magenta} \text{stable}}
  \\
  (iv) & \|R_w(\lambda) \| \approx 0 \;\; \text{with} \;\; R_w(\lambda) \;\; {\color{magenta} \text{stable}}
  \end{array}
```
When solving an AFDIP, it is usually attempted to fulfil the stronger requirement that $\overline R_f(\lambda) =  0$ (which is equivalent to require $S_{R_f} = S$). If this is not feasible, then the above problem is solved.

For a partitioned filter corresponding to a targeted structure matrix $S$, we can define for the $i$-th filter component the associated value of the fault-to-noise gap, which  characterizes the noise attenuation properties of the $i$-th filter. Let $f^{(i)}$ be formed from the subset of faults corresponding to nonzero entries in the $i$-th row of $S$ and  let $\bar f^{(i)}$ be formed from the complementary subset of faults corresponding to zero entries in the $i$-th row of $S$. If $R_{f^{(i)}}^{(i)}(\lambda)$ and $R_{\bar f^{(i)}}^{(i)}(\lambda)$  are formed from the columns of $R_{f}^{(i)}(\lambda)$ corresponding to $f^{(i)}$ and $\bar f^{(i)}$, respectively, then the fault-to-noise gap of the $i$-th filter can be defined as
```math 
 J_{2}^{(i)} =   \big\| R_{f^{(i)}}^{(i)}(\lambda) \big\|_{\infty -} / \big\|\big[\, R_{\bar f^{(i)}}^{(i)}(\lambda) \;  R_{w}^{(i)}(\lambda)\,\big]\big\|_\infty . 
```
For a similar characterization of the strong complete fault detectability of the subset of faults corresponding to nonzero entries in the $i$-th row of $S$,  we have
```math 
 \widetilde J_{2}^{(i)} =   \big\| R_{f^{(i)}}^{(i)}(\lambda) \big\|_{\Omega -} / \big\|\big[\, R_{\bar f^{(i)}}^{(i)}(\lambda) \;  R_{w}^{(i)}(\lambda)\,\big]\big\|_\infty . 
```

For the synthesis of fault detection and isolation filters which solve the AFDIP the function [`afdisyn`](@ref) is available and for the evaluation of the  fault sensitivity condition the function [`fdif2ngap`](@ref) is available. 


### Exact model-matching problem (EMMP)

Let $M_r(\lambda)$ be a given $q\times m_f$ TFM of a stable reference model  specifying the desired input-output behavior from the faults to residuals as
```math 
{\mathbf{r}}(\lambda) = M_r(\lambda) {\mathbf{f}}(\lambda). 
```
Together with the decoupling conditions $R_u(\lambda) = 0$ and $R_d(\lambda) = 0$, the determination of  $Q(\lambda)$ involves the solution of a linear matrix equation with rational function coefficients. However, a particular choice of $M_r(\lambda)$ may lead to a solution $Q(\lambda)$ which is not proper or is unstable or has both these undesirable properties. Therefore, besides determining $Q(\lambda)$, 
the determination of a suitable updating factor $M(\lambda)$ of $M_r(\lambda)$ is necessary to ensure the stability of the solution $Q(\lambda)$ for $R_f(\lambda) = M(\lambda) M_r(\lambda)$ (and also of $R_w(\lambda)$). This leads to the following formulation of the EMMP:

**EMMP:** Determine a stable residual generator $Q(\lambda)$ and a 
stable, diagonal, and invertible $M(\lambda)$ such that 
```math 
\begin{array}{rl} (i) & R_u(\lambda) = 0 \\ (ii) & R_d(\lambda) = 0 \\
  (iii) & R_f(\lambda) = M(\lambda)M_r(\lambda) \\
  (iv) &R_w(\lambda) \;\; {\color{magenta} \text{stable}} 
  \end{array}
```

A typical choice for $M_r(\lambda)$ is an $m_f \times m_f$  diagonal and invertible TFM, which ensures that each residual $r_i(t)$ is influenced only by the fault $f_i(t)$. This would allow the isolation of arbitrary combinations of up to $m_f$ simultaneous faults.  The choice $M_r(\lambda) = I_{m_f}$ targets the solution of an  _exact fault estimation problem_ (EFEP).


For the synthesis of fault detection and isolation filters which solve the EMMP the function [`emmsyn`](@ref) is available. This function can also address the solution of the EMMP with more general reference models (e.g., having components from the control inputs and/or disturbance inputs). 

### Approximate model-matching problem (AMMP}

Let $M_r(\lambda)$ be a given $q\times m_f$ TFM of a stable reference model  specifying the desired input-output behavior from the faults to residuals as
```math 
{\mathbf{r}}(\lambda) = M_r(\lambda) {\mathbf{f}}(\lambda). 
```
**AMMP:** Determine a stable residual generator $Q(\lambda)$ and a 
stable, diagonal, and invertible $M(\lambda)$ such that 
```math 
\begin{array}{rl} (i) & R_u(\lambda) = 0 \\ (ii) & R_d(\lambda) = 0 \\
  (iii) & R_f(\lambda) \approx M(\lambda)M_r(\lambda) \\
  (iv) & \|R_w(\lambda) \| \approx 0 \;\; \text{with} \;\; R_w(\lambda) \;\; {\color{magenta} \text{stable}}
  \end{array}
```

A criterion suitable to characterize the solution of approximate model-matching based syntheses is the  residual error norm
```math 
J_3 = \big\| R(\lambda)- M_r(\lambda)\big\|_{\infty/2}, 
```
where $R(\lambda) = R_f(\lambda)$ and $M_r(\lambda)$ the reference model (possibly updated). For more generality, this criterion can be defined with 
$R(\lambda) = [\, R_u(\lambda)\; R_d(\lambda)\; R_f(\lambda)\; R_w(\lambda) \,\,]$, the resulting internal form, and $M_r(\lambda)$ the desired reference model $M_r(\lambda) = [\, M_{ru}(\lambda)\; M_{rd}(\lambda)\; M_{rf}(\lambda)\; M_{rw}(\lambda)\,]$. When applied to the results computed by other synthesis approaches (e.g., to solve the AFDP or AFDIP), this criterion can be formulated as
```math 
\widetilde J_3 = \big\| R_w(\lambda)\big\|_{\infty/2}, 
```
which corresponds to assume that $M(\lambda) = I$ and $M_r(\lambda) = [\, R_u(\lambda)\; R_d(\lambda)\; R_f(\lambda)\; 0 \,]$ (i.e., a perfect matching of control, disturbance and fault channels is always achieved).

For the synthesis of fault detection and isolation filters which solve the AMMP the function [`ammsyn`](@ref) is available. This function can also address the solution of the AMMP with more general reference models (e.g., having components from the control inputs, disturbance inputs and noise inputs). For the evaluation of the  model-matching performace the function [`fdimmperf`](@ref) is available. 




### References

[1]   A. Varga, Solving Fault Diagnosis Problems – Linear Synthesis Techniques, Vol. 84 of Studies in Systems, Decision and Control, Springer International Publishing, 2017.



