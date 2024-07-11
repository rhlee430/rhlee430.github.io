---
title: Exact Diagonalization of Composite Fermi Liquids
---
In this note, we deal with exact diagonalization of composite fermi liquids. We first construct Bloch state bases of Landau levels, and subsequently write the Hamiltonian and diagonalize it in such bases.

## Bloch Wave Functions

We follow Ref. [^1] for preparing reference states for the diagonalization. It is a way that Ref. [^2] conducted ED. First we construct torus spanned by to vectors $\textbf{L}_1= N_1 \textbf{a}_1, \textbf{L}_2= N_2 \textbf{a}_2$. Note that the torus is yet continuum. Let corresponding reciprocal lattice vectors $\textbf{b}_1$ and $\textbf{b}_2$. Assume that total number of flux is  $N_\phi=N_1N_2$, i.e. a single flux is coupled for each unit cell. Starting from a strip basis $|j\rangle$, we could construct new bases characterized by conserved pseudo momentum $\textbf{k}=k_1 \textbf b_1+k_2 \textbf b_2 \in BZ$,

$$
|\textbf{k}\rangle = |k_1, k_2\rangle = \frac{1}{\sqrt{N_1}} \sum_{m=0}^{N_1-1} e^{i2\pi m k_1/N_1} |j=m N_2 + k_2 \rangle.
$$

These states are periodic in $k_1$, $|k_1+N_1, k_2\rangle =|k_1, k_2\rangle$, but *quasi*-periodic in $k_2$, $|k_1, k_2+N_2\rangle = e^{-i2\pi k_1/N_1} |k_1, k_2\rangle$. The LLL-projected density is 
$$
\rho_{\textbf{q}}=e^{-\textbf{q}^2 l^2/4} \sum_{\textbf{k}}^{BZ} e^{-i 2\pi q_1(k_2+q_2/2)/N_\phi} |\textbf{k}\rangle \langle \textbf{k}+\textbf{q}|
$$
LLL-projected Interactions are able to be written in these bases,
$$
\begin{align}
H_{int} &=\sum_{\textbf{q}}V_{\textbf{q}} :\rho_{\textbf{q}}\rho_{-\textbf{q}}: \\
&=\sum_{\textbf{k}, \textbf{k}^\prime, \textbf{q}} e^{i 2\pi q_1(k_2-k^\prime_2+q_2)/N_\phi} V_{\textbf{q}} e^{-\textbf{q}^2 l^2/2} c^\dagger_{\textbf{k}+\textbf{q}} c^\dagger_{\textbf{k}^\prime-\textbf{q}} c_{\textbf{k}^\prime} c_{\textbf{k}}.
\end{align}
$$
As we are dealing with [[the Haldane pseudopotential]] of nLL, we should substitute $v^{(n)}_{\textbf{q}}$ in $V_\textbf{q}$.

---
[^1]: Y.-L. Wu, N. Regnault, and B. A. Bernevig, _Bloch Model Wave Functions and Pseudopotentials for All Fractional Chern Insulators_, Phys. Rev. Lett. **110**, 106802 (2013).
[^2]: H. Goldman, A. P. Reddy, N. Paul, and L. Fu, _Zero-Field Composite Fermi Liquid in Twisted Semiconductor Bilayers_, Phys. Rev. Lett. **131**, 136501 (2023).

