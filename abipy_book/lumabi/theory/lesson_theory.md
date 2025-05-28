---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.10.3
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Theory

This chapter presents the theoretical formalism behind the computation of the luminescence properties of point defects in solids.
It starts with a section dedicated to the theory of luminescence, showing how Fermi's golden rule is approximated
to obtain a tractable way of computing the luminescence lineshape, following either an effective phonon mode model or a multi-phonon mode model.
Then, the computational methodology used to obtain the parameters entering the lineshape equation is presented.
Subtleties regarding the increase of the phonon supercell size thanks to the use of the forces instead of displacements are discussed.
Finally, the IFCs embedding approach, that allows to compute the phonon modes of a defect system in large supercells, is presented.

## From the Fermi's Golden rule to the luminescence lineshape

Following Fermi's golden rule, the absolute luminescence intensity $I(\hbar\omega)$ (number of photons per unit time per unit energy)
associated to one emitting center with two states e and g, is expressed as a function of the photon energy $\hbar\omega$ as

$$
    I(\hbar \omega)=\frac{n_D \omega^3}{3 \pi \hbar \epsilon_0 c^3} |\boldsymbol{\mu_{eg}}|^2 \delta\left(E_g-E_e-\hbar \omega\right),
$$ (first_lum_intensity)

where $n_D$ is the material refractive index,
$\boldsymbol{\mu_{eg}}=\langle\Psi_{e}|\boldsymbol{\mu}|\Psi_{g}\rangle$ is the total dipole matrix element.
It is then assumed that the electronic part of the dipole matrix element between the excited and ground state depends
weakly on the nuclear coordinates (Franck-Condon approximation).
We also suppose that the nuclear motions can be expressed as a superposition of 3N harmonic normal modes of vibration $\nu$,
represented by normal coordinates $Q_{\nu}$ and frequency $\omega_{\nu}$, with N being the number of atoms in the cell considered.
The total vibrational state $\chi_{\boldsymbol{n}}$ is the product of 3N independent harmonic
oscillator eigenfunctions $\chi_{\boldsymbol{n_{\nu}}}$ with $n_{\nu}$ the vibrational state of the $\nu$-th harmonic oscillator.
We write the normalized luminescence intensity as

$$
L(\hbar\omega)=C\omega^3A(\hbar\omega),
$$ (norm_lum_intensity)

with C a normalization constant.
In the case of the absolute intensity $I(\hbar\omega)=C\omega^3A(\hbar\omega)$,
$C=\frac{n_D}{3 \pi \hbar \epsilon_0 c^3}|\boldsymbol{\mu}_{eg}^{el}|^2$
with $\boldsymbol{\mu}_{eg}^{el}$ the purely electronic dipole moment (inversely proportional to the radiative lifetime of the emitting center).
The emission lineshape function $ A(\hbar \omega)$, which dictates the shape of the spectrum, writes

$$
    A(\hbar \omega)= \sum_{\boldsymbol{n,m}}p_{m_{\nu}}(T)|\langle\chi_{e,\boldsymbol{m}}|\chi_{g,\boldsymbol{n}}\rangle|^2 \delta\left(E_{ZPL}-E_{g,\boldsymbol{n}}+E_{e,\boldsymbol{m}}-\hbar \omega\right).
$$ (second_lum_intensity)

In equation {eq}`second_lum_intensity`, $\boldsymbol{n}$ denotes the set of 3N vibrational states {$n_1$,$n_2$,...,$n_{3N}$}, $E_{ZPL}$
is the so-called zero-phonon line energy (purely electronic energy difference),
$E_{g/e,\boldsymbol{n/m}}=\sum{_\nu}n_{\nu}\hbar\omega_{g/e,\nu}$ is the vibrational energy of state $\chi_{g/e,\boldsymbol{n/m}}$,
and $p_{_{\nu}}(T)$ is a Bose-Einsten thermal occupation factor.

The vibrational modes of the excited and ground states can be in principle different, making equation {eq}`second_lum_intensity`
computationally very expensive because of the highly multidimensional integrals to be evaluated (Franck-condon overlaps).
Some approximations should be made.
We assume first that the vibrational modes in the excited electronic state are identical to those in the ground state.
We also assume that we work at T=0 K, to alleviate the equations ($p_{m_{\nu}}(T=0)$ = 0 if $m_{\nu}$ a vibrational state
of the electronic excited state is different from 0).
Considering temperature dependent luminescence spectrum will be presented later.

### Effective phonon mode model

The simplest model to start with is built when considering that the 3N phonon modes of the system can be reduced to a single effective phonon mode.
Suppose that a single effective mode couples to the electronic transition,
the Franck-Condon overlaps will reduce to $|\langle\chi_{e,0}|\chi_{g,n}\rangle|^2$.
In the harmonic approximation, the vibrational eigenfunctions $\chi$ are expressed with Hermite polynomials and
their overlaps can be computed analytically :

$$
|\langle\chi_{e,0}|\chi_{g,n}\rangle|^2=e^{-S}\frac{S^n}{n!},
$$ (HuangRhys_eq)

with $S$ the Huang-Rhys factor, a dimensionless constant giving the mean number of effective phonons of frequency $\Omega$
involved in the electronic transition,

$$
S=\frac{\frac12\Omega^2\Delta Q^2}{\hbar\omega}=\frac{\Omega\Delta Q^2}{2\hbar},
$$ (HuangRhys_factor)

and $\Delta Q$ the offset between the harmonic potential energy surface of the ground and excited state.
In practice, such offset is computed as a total mass weighted atomic displacements induced by the electronic transition:

$$
(\Delta Q)^2=\sum_{\kappa i}m_\kappa(R_{\text{e};\kappa i}-R_{\text{g};\kappa i})^2,
$$ (Delta_Q)

where $i$ labels the cartesian axes, $\kappa$ the atoms, $m_{\kappa}$ atomic masses, $R_{e;\kappa i}$
and $R_{g;\kappa i}$ are respectively atomic positions in the excited and the ground states.

This allows to simplify equation {eq}`second_lum_intensity` to

$$
A(\hbar\omega)=\sum_n e^{-S}\frac{S^n}{n!}\delta(E_{\mathrm{ZPL}}-n\hbar\Omega-\hbar\omega).
$$ (lineshape_1D)

A schematic figure of such model is presented on figure {numref}`displaced_HA` for a single mode $\nu$.

```{figure} displaced_HA.png
---
height: 300px
name: displaced_HA
---
Schematic representation of the origin of photoluminescence (PL) spectra.
On the left, ground state (GS) and excited state (ES) energy curves are projected along phonon normal coordinate $Q_{\nu}$ (here only for one mode)
and are approximated by harmonic functions with the same frequency $\omega_{\nu}$.
Vibrational energy levels and corresponding eigenfunctions are shown as horizontal lines and colored areas.
GS and ES are displaced by $\Delta Q_{\nu}$, the mass-weighted displacement between the minimum of the GS and the ES curves.
On the right, the PL spectrum is formed.
The zero-phonon line (ZPL) comes from the transition between the first vibrational level of the ES to the first vibrational level of the GS.
Other transitions give the phonon sideband (PSB).
The intensity of each peak is computed with the overlap between corresponding eigenfunctions.
```

In order to visualise the effect of harmonic oscillator frequency or the offset $\Delta Q$ on the luminescence lineshape,
you can slide the cursors below.
For the frequency, we assume harmonic oscillators with frequency $\omega=\sqrt{\frac{k}{\mu}}$ with $k$ a force/spring constant,
and $\mu$ the oscillator mass.


```{code-cell}
from lineshape_slider import get_plotly_lineshape
%matplotlib inline
```

```{code-cell}
get_plotly_lineshape(slider="frequency", num_step=10)
```

```{code-cell}
get_plotly_lineshape(slider="Delta_Q", num_step=10)
```

### Multi phonon mode model and generating function

We now consider the 3N vibrational modes.
Quantities of the above section are now computed "per phonon mode".
The $\Delta Q_{\nu}$ (mass-weighted atomic displacement projected along phonon mode $\nu$)

$$
\Delta Q_\nu=\sum_{\kappa\alpha}\sqrt{m_\kappa}\Delta R_{\kappa\alpha}e_{\nu,\kappa\alpha}
$$ (Delta_Q_multi)

are used to compute the partial Huang-Rhys factor $S_\nu=\frac{\omega_\nu\Delta Q_\nu^2}{2\hbar}$.
$e_{\nu,\kappa\alpha}$ are the phonon eigenvectors.
These define the Huang-Rhys spectral function:

$$
S(\hbar\omega) = \sum_{\nu}S_{\nu}\delta(\hbar\omega-\hbar\omega_{\nu}).
$$ (HuangRhys_spectral_function)

When dealing with 3N phonons, a direct evaluation of equation {eq}`second_lum_intensity` is impractical.
Instead, the so-called generating function approach is used.
The lineshape function $A(\hbar\omega)$ is evaluated as the Fourier transform of the generating function $G(t)$

$$
A(\hbar\omega,T)=\int_{-\infty}^{+\infty}G(t,T)e^{i\omega t-\frac\gamma\hbar|t|-i\frac{E^\mathrm{ZPL}}\hbar t}dt,
$$ (lineshape_generating)
with

$$
G(t,T=0)=e^{S(t)-S(0)},
$$ (generating_fct)

where $S(t)=\sum_\nu S_\nu e^{i\omega_\nu t}$ is the Fourier transform of the Huang-Rhys spectral function.
$\gamma$ is a constant homogeneous Lorentzian broadening associated to each vibronic transition.

While working in the time domain, the effect of temperature (transitions involving initial vibrational states $n_{\nu} \ne 0$,
weighted by Bose-Einstein occupation) can be included straightforwardly by rewriting the generating function

$$
G(t,T)=e^{S(t)-S(0)+C(t,T)+C(-t,T)-2C(0,T)}
$$ (generating_fct_T)

with $C(t,T)=\sum_\nu\overline{n}_\nu(T)S_\nu e^{i\omega_\nu t}$ the Fourier transform of the temperature weighted Huang-Rhys spectral function,
and $\overline{n}_\nu(T)$ the average occupation number of $\nu$-th phonon mode:

$$
\overline{n}_\nu(T)=\frac{1}{e^{\frac{\hbar\omega_\nu}{k_BT}}-1}.
$$ (bose_einstein)

One can connect the multi-phonon modes methodology to the simpler effective phonon mode model presented in the previous section.
The total normal coordinate change $\Delta Q$, due to the orthonormality of the phonon eigenvectors is linked through
the partial $\Delta Q_{\nu}$ through:

$$
(\Delta Q)^2=\sum_{\nu}(\Delta Q_{\nu})^2,
$$ (Delta_Q_total)

and allows one to define the weight by which a mode contributes to the total atomic relaxation:

$$
p_{\nu}=(\Delta Q_{\nu}/\Delta Q)^2.
$$ (p_nu)

It is then possible to define an effective frequency as

$$
\omega_{\mathrm{eff}}^2=\sum_{\nu}p_{\nu}\omega_{\nu}^2,
$$ (omega_eff)

and the total Huang-Rhys factor as:

$$
S=\sum_{\nu}S_{\nu}.
$$ (HuangRhys_total)

### Semi-classical approach

Within a semi-classical formulation {cite}`henderson2006optical`, one can find formulas for the full width a half maximum of the emission shape:

$$
W(0)=S_{\mathrm{em}}\hbar\Omega_{\mathrm{g}}\sqrt{8\ln2}/\sqrt{S_{\mathrm{abs}}}.
$$ (fwhm_semi_classical)

where $S_{\mathrm{em}}$ and $S_{\mathrm{abs}}$ are the Huang-Rhys factors associated to the emission and absorption processes respectively,
and $\Omega_{\mathrm{g}}$ is the effective frequency of the ground state.

Averaging the Bose-Einstein occupations of the initial vibrational states allows one to compute the temperature dependent FWHM:

$$
W(T)=W(0)\sqrt{\coth(\hbar\Omega_\mathrm{e}/2k_BT)},
$$ (fwhm_semi_classical_T)

## Computational methodology

From the above approximations, we see that the photo-luminescent lineshape of a point defect in a solid is uniquely
defined by the Huang-Rhys spectral function $S(\hbar\omega)$ and the zero-phonon line.
In order to compute it from first-principles, we need to obtain:

-  The zero-phonon line energy $E_{ZPL}$, which is the total energy difference between the electronic excited and ground states of the system.
-  The atomic relaxation associated to the change in electronic state, obtained as the difference between the relaxed atomic positions
   of the excited state and relaxed atomic positions of the ground state.
-  The vibrational modes of the system. Note that in the case of the effective phonon model,
   this is not required since the effective mode is obtained from the atomic relaxation.

One way to obtain the first two parameters is to use DFT following the $\Delta SCF$ constrained occupation method.
The vibrational modes can be obtained with finite difference or DFPT, or by using the embedding methodology presented in latter in this chapter.


### $\Delta$SCF constrained occupation method

This method refers to the use of DFT with non-Aufbau electronic occupations to mimick the electron-hole interaction.
Transition energies are computed by taking the differences between total DFT energies with different occupations.
Note that the excited state occupations are specific to the system under study.
In the case of Eu$^{2+}$, one of the seven 4f electron of the spin-up channel is promoted to the next spin-up 5d energy state,
as shown in {numref}`delta_scf`.
This figure also illustrates the connection with the configuration coordinate model.


```{figure} delta_scf.png
---
width: 100%
name: delta_scf
---
Schematic representation of the constrained occupation $\Delta$SCF method, as used in this work for simulating the luminescent properties of Eu$^{2+}$ phosphor.
```

The procedure is as follows:

- Start with a relaxed ground state (labeled $A_g$) with energy $E_{g}$.
- Promote an electron to obtain the excited state configuration ($A_g^*$, energy $E_{g}^{*}$).
  The difference in energy with the $A_g$ state provides an estimate of the semi-classical absorption energy:

  $$
  E_{\mathrm{abs}} = E_{g}^{*} - E_{g}
  $$ (absorption_energy)

- After structural relaxation in the excited state, obtain the relaxed excited state ($A_e^*$, energy $E_e^*$), which allows computation of the zero-phonon line (ZPL) energy:

  $$
  E_{\mathrm{ZPL}} = E_e^* - E_g
  $$ (ZPL_energy)
- The Franck-Condon relaxation energy of the excited state is:

  $$
  E_{\mathrm{FC,e}} = E_{g}^{*} - E_{e}^{*}
  $$ (FC_e)

- De-promoting the electron with the relaxed atomic positions of the excited state yields the ground state with excited geometry ($A_e$, energy $E_e$). The emission energy is:

  $$
  E_{\mathrm{em}} = E_e^* - E_e
  $$ (emission_energy)

- The Franck-Condon relaxation energy of the ground state is:

  $$
  E_{\mathrm{FC,g}} = E_{e} - E_{g}
  $$ (FC_g)

Assuming harmonicity, the effective frequencies of the ground or excited states can be extracted as

$$
\omega_{\mathrm{eff},\{g,e\}}^2 = \frac{2 E_{\mathrm{FC},\{g,e\}}}{\Delta Q^2}
$$ (omega_eff_g_e)

where $\Delta Q$ is the mass-weighted atomic displacement between the ground and excited states, as defined in equation {eq}`Delta_Q`.
The corresponding Huang-Rhys factors can be obtained from {eq}`HuangRhys_factor`.

### Forces vs Displacements, increase of supercell size

For a given phonon mode $\nu$, the corresponding partial Huang-Rhys factor $S_{\nu}=\frac{\omega_{\nu}\Delta Q_{\nu}^2}{2\hbar}$
is computed thanks to the mass-weighted displacement between ground and excited states projected on this phonon mode:

$$
	\Delta Q_\nu=\sum_{\kappa}\sqrt{M_{\kappa}}\Delta R_{\kappa}{e_{\nu,\kappa}},
$$ (Delta_Q_nu_emb)

where $\kappa$ labels the atoms and $\alpha$ the cartesian direction.
Under the harmonic approximation,

$$
	M_{\kappa}\omega_{\nu}^2e_{\nu,\kappa}=\sum_{\kappa'}C_{\kappa,\kappa'}e_{\nu,\kappa'},
$$ (IFCs)

with the interatomic force constants (IFCs) $C_{\kappa,\kappa'}=\frac{\Delta F_{\kappa'}}{\Delta R_{\kappa}}$, Eq. {eq}`Delta_Q_nu_emb` rewrites:

$$
	\Delta Q_\nu=\frac{1}{\omega_\nu^2}\sum_{\kappa}\frac{\Delta F_{\kappa} e_{\nu,\kappa}}{\sqrt{M_{\kappa}}}.
$$ (Delta_Q_nu_emb_forces)

The use of Eq.{eq}`Delta_Q_nu_emb_forces` allows increasing the supercell size and hence minimize finite-size effects.

Indeed, the forces decay faster to zero compared to the displacements with respect to the distance from the defect
as illustrated in see {numref}`Forces_vs_dis`, panel **a-b**.

We assume that the forces computed within the small red supercell (here 288 atoms), see {numref}`Forces_vs_dis`, panel **c** ,
are already converged and without finite-size effect because of this rapid decay, and that the forces
outside this red supercell are essentially zero.

This means that we can use the forces computed in this red supercell, and deduce the displacements at much larger distance
in the blue supercell via the IFCs computed within this large blue supercell.

This subtle point is actually directly included when using Eq.{eq}`Delta_Q_nu_emb_forces` in the blue supercell.
This also means that the IFCs/phonons modes should be computed in the same blue supercell.
This is discussed in the next section.

```{figure} Forces_vs_dis.png
---
height: 300px
name: Forces_vs_dis
---
(a) Norm of the displacement induced by an electronic transition as a function of the distance from the defect atom (Eu).
(b) Norm of the forces in the ground state at the excited state equilibrium atomic positions. The decay of forces (b)
is much faster than the decay of displacements (a).
(c) Cartoon of a small red supercell containing substitutional defect where the forces are the ones computed with DFT.
Outside this red supercell, the forces are set to zero because of their short-range decay while the displacements
are non-zero and are computed with the IFCs computed in this large blue supercell.
```

### IFCs embedding

In the computation of the partial Huang-Rhys factors, one would like to obtain the phonons in very large supercells
at wave vector $\mathbf{q}\rightarrow(0,0,0)$ ($\Gamma$ point), so that the coupling with long-wavelength phonons is correctly captured.
One would also like to include the coupling with (localized) phonons modes introduced by the defect.
A direct approach, either with DFPT or finite difference, is not computationally attractive.
Indeed, DFPT is best used for small pristine primitive cells at dense $\mathbf{q}$-mesh
(which corresponds to large pristine supercells after a folding procedure), while finite difference, for which the introduction
of defect does not cause problem, remains costly for large supercells.
One way to include both the defect effect (local modes) while converging long-wavelength phonons is to employ
the IFCs embedding approach {cite}`alkauskas2014,jin2021photoluminescence`.


```{figure} emb_approach.png
---
height: 300px
name: emb_approach
---
Schematic view of the interatomic force constants (IFC) embedding approach to obtain defect phonons in large supercells.
Each arrow represents an IFC between a pair of two atoms. See text for details.
```

The procedure, illustrated in figure {numref}`emb_approach`, is as follows:

First, the real-space interatomic force constants (IFCs) of the defect system, denoted as $C_{\kappa\alpha,\kappa’\alpha}^{\mathrm{defect}}$,
are computed within a relatively small supercell (typically with a few hundred atoms).
This is often done using a finite difference approach on the supercell at the $\Gamma$ point,
though Density Functional Perturbation Theory (DFPT) could also be employed if the memory requirements are manageable and/or
if the level of theory required to describe the defect (e.g., DFT+U) is supported by the DFPT code.

Second, the IFCs of the pristine system, $C_{\kappa\alpha,\kappa’\alpha}^{\mathrm{pristine}}$, are obtained within
a much larger supercell (typically with a few thousand atoms).
A practical approach is to use DFPT on the unit cell to calculate the dynamical matrix on a $\mathbf{q}$-mesh,
potentially followed by Fourier interpolation to obtain the matrix on a denser [$N_x,N_y,N_z$] $\mathbf{q}$-mesh.
The information in the dynamical matrix, computed on a unit cell for a dense [$N_x,N_y,N_z$] $\mathbf{q}$-mesh,
is equivalent to that in a dynamical matrix for a [$N_x,N_y,N_z$] supercell at $\mathbf{q} = [0,0,0]$ (the $\Gamma$ point).
Mapping a $\mathbf{q}$-mesh to $\mathbf{q} = [0,0,0]$ is referred to a folding procedure,
and yields the pristine IFCs $C_{\kappa\alpha,\kappa’\alpha}^{\mathrm{pristine}}$.

Third, an embedded IFC matrix $C_{\kappa\alpha,\kappa’\alpha}^{\mathrm{emb}}$ is constructed using the following rules:
If both atoms $\kappa$ and $\kappa’$ are within a sphere centered around the defect, with a cut-off radius $R_c$,
then $C_{\kappa\alpha,\kappa’\alpha}^{\mathrm{emb}}$ is set to $C_{\kappa\alpha,\kappa’\alpha}^{\mathrm{defect}}$.
For all other atomic pairs, the embedded IFC is set to the pristine value: $C_{\kappa\alpha,\kappa’\alpha}^{\mathrm{emb}} = C_{\kappa\alpha,\kappa’\alpha}^{\mathrm{pristine}}$.
The cut-off radius $R_c$ is typically determined by the size of the initial supercell used to compute the defect IFCs.
Additionally, one could use a second cut-off radii $R_b$.
If two atoms are not within a sphere of radii $R_b$, the IFC is set to zero.
The use of this second cut-off has been used in some studies {cite}`razinkovas2021vibrational,jin2021photoluminescence`
in order to obtain large sparse matrices, for which dedicated techniques allows faster diagonalization,
see for instance reference {cite}`razinkovas2021vibrational`.
In this work, we do not use this second cut-off.

Overall, the whole procedure can break the acoustic sum rule

$$
C_{\kappa\alpha,\kappa\beta} = - \sum_{\kappa'\ne\kappa}C_{\kappa\alpha,\kappa'\beta}.
$$ (acoustic_sum_rule)

To correct for this, we set:

$$
C_{\kappa\alpha,\kappa\alpha}^{\mathrm{emb}}=-\sum_{\kappa'\ne\kappa}C_{\kappa\alpha,\kappa'\alpha}^{\mathrm{emb}}.
$$ (acoustic_sum_rule_emb)

For the treatment of the Born effective charges (BECs), if the BEC of the substituted atom was computed in the defect phonon calculation,
it is replaced by the BEC of the dopant atom.
If not, one can choose to use the most common oxidation state of the dopant atom as its new BEC.
The same approach applies to interstitials, but in this case, the dopant atom and its corresponding BEC are simply added.
For vacancies, the atom and its BEC are removed.

Diagonalizing this embedded IFCs provides the phonon modes of a defectuous [$N_x,N_y,N_z$] supercell, where the defect effect on the IFCs
was included through the embedding procedure.

The localization of the phonon modes can be characterized by computing the inverse participation ratio (IPR)
defined as {cite}`alkauskas2014,jin2021photoluminescence`

$$
\mathrm{IPR}_{\nu}=\frac{1}{\sum_{\kappa}|\langle \mathbf{e}_{\nu,\kappa}| \mathbf{e}_{\nu,\kappa} \rangle|^2},
$$ (IPR)

where $\mathbf{e}_{\nu,\kappa}$ is the eigenvector of the phonon mode $\nu$ associated to atom $\kappa$.
The IPR is a measure of the localization of the phonon mode $\nu$ and is normalized to the number of atoms in the supercell $N$.
The IPR can be interpreted as a measure of the number of atoms that participate in a given phonon mode,
which indicates, roughly speaking, the number of atoms that participate to a phonon mode $\nu$.
For example, $\mathrm{IPR}_{\nu}=1$ means that only one atom vibrates and therefore the phonon mode is maximally localized,
while $\mathrm{IPR}_{\nu}=N$ means that N atoms vibrates in the supercell with the same amplitude.

The localization ratio $\beta_{\nu}$ is defined by taking the ratio of the total number of atoms
in the supercell $N$ and the $\rm{IPR}$, $\beta_{\nu} =N/\mathrm{IPR}_{\nu}$
where $\beta_{\nu} \approx 1$ represents a bulk-like delocalized mode while $\beta_{\nu} \gg 1$ corresponds to a quasi-local or local mode.


## References
```{bibliography}
:style: unsrt
```
