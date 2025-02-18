A wave guide model for a transverse magnetic (TM) solution

The case is defined such that the time-harmonic solution is expected to be of
the form

 $$ \mathbf{E} = (\mathrm{cos}(k_1 x), 0, -ik_1/k_3 \mathrm{sin}(k_1 x)) E_0 \mathrm{exp}[-i(w t - k_3 z)].  $$

This case is checked to give the relative ($\mathbf{L}_2$) error

$$ || \mathbf{E} - \mathbf{E}_h || / || \mathbf{E} || = 4.62445 \cdot 10^{-2}. $$

If we apply a regular subdivision (Mesh Levels = 2), we obtain

$$ || \mathbf{E} - \mathbf{E}_h || / || \mathbf{E} || = 1.17795 \cdot 10^{-2} $$

so we see quadratic convergence as expected.