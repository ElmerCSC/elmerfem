A wave guide model for a transverse magnetic (TM) solution

The case is defined such that the time-harmonic solution is expected to be of
the form

 $$ \mathbf{E} = (\mathrm{cos}(k_1 x), 0, -ik_1/k_3 \mathrm{sin}(k_1 x)) E_0 \mathrm{exp}[-i(w t - k_3 z)].  $$

This case has been used to verify that the convergence rates which are expected for different finite element approximations are obtained in practice. Here Nd_1(p=k) and Nd_2(p=k) refer to the Nedelec approximation of the first and second kind, respectively, and k is the order of basis functions. 


|     | Nd_1(p=1)      | Nd_2(p=1)      | Nd_1(p=2)      | Nd_2(p=2)      |
|-----|----------------|----------------|----------------|----------------|
| h   | 0.11247701E+00 | 0.45481827E-01 | 0.56386830E-02 | 0.48199766E-03 |
| h/2 | 0.56636760E-01 | 0.11677471E-01 | 0.14273238E-02 | 0.58223424E-04 |
|     | rate r = 1     | r = 2          | r = 2          | r = 3          |

**Table I:** The relative ($\mathbf{L}_2$) error $|| \mathbf{E} - \mathbf{E}_h || / || \mathbf{E} ||$ for
different finite element approximations.


|     | Nd_1(p=1)      | Nd_2(p=1)      | Nd_1(p=2)      | Nd_2(p=2)      |
|-----|----------------|----------------|----------------|----------------|
| h   | 0.10852124E+00 | 0.11647564E+00 | 0.76389743E-02 | 0.76299899E-02 |
| h/2 | 0.54456111E-01 | 0.55530316E-01 | 0.19122851E-02 | 0.19115754E-02 |
|     | rate r = 1     | r = 1          | r = 2          | r = 2          |

**Table II:** The relative ($\mathbf{L}_2$) error $|| \mathbf{curl}(\mathbf{E} - \mathbf{E}_h) || / || \mathbf{curl}\mathbf{E} ||$ for different finite element approximations.
