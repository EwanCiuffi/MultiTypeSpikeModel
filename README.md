# Multi-type Spike Model

The **Multi-type Spike Model** is a BEAST 2 clock model for modelling **punctuated evolution** under structured birth–death processes.

It extends the Gamma Spike model of Douglas et al. (2025) by allowing:

- **State‑dependent speciation and extinction** 
- **Time varying (skyline) phylodynamic parameters**

> **Note:** Unlike the original Gamma Spike model, the Multi-type Spike Model **does not have the option to estimate the number or timing of hidden events** (“stubs”) explicitly. It integrates over the number of hidden speciation events numerically. However, all other features of the Gamma Spike model are retained.

---

## Model overview

The evolutionary distance along branch $e$ is:

$$
d^e = r^e \ \tau^e + s^e
$$

where:

- $r^e$ is the **branch rate**  
- $\tau^e$ is the **branch length** (time)  
- $s^e$ is the **spike** on branch $e$
  
The term $r^e \tau^e$ represents gradual change along the branch, while $s^e$ represents the contribution of punctuated evolution (“spikes”) along the branch.

---

## Spike Distribution

Each branching event contributes an instantaneous burst of substitutions.
If branch $e$ belongs to type $i$ and has $N_i^e$ branching events
(observed + hidden), then its total spike magnitude is:

$$
s_i^e \sim S^\mu_i \cdot \mathrm{Gamma}\\left(
\text{shape} = S^\alpha \cdot N_i^e\ , 
\text{scale} = \frac{1}{S^\alpha}
\right).
$$

---

## Model Averaging
As in the original Gamma Spike model, the posterior support for punctuated evolution can be assessed directly from the data. The clock indicator (`useSpikeModel`) is assigned the prior

$\mathbb{I} ∼ Bernoulli(0.5)$,

such that the posterior probability of $\mathbb{I} = 1$ provides evidence for punctuated evolution relative to a relaxed (or strict) clock model.

---

### Key Parameters

| Parameter | Meaning | Notes |
|----------|---------|-------|
| $S^\mu$ or $S^\mu_i$ | Spike mean | Expected substitutions per site **per branching event**. Type‑specific in multi‑type mode. |
| $S^\alpha$ | Spike shape | Controls variance in spike sizes. Larger $S^\alpha$ $\rightarrow$ more uniform spikes. |
| $\mathbb{I}$ | Spike model indicator | Bernoulli prior; posterior support quantifies evidence for punctuated evolution. |


> In single‑type models, all lineages share the same $S^\mu$.  
> In multi‑type models, each type $i$ can have its own $S^\mu_i$, enabling type‑specific punctuated evolution.

---

## Testing for Differences in Punctuated Evolution Between Types

Type‑specific spike means allow for hypothesis testing:

- **Trait‑dependent evolution:**  
  Are certain traits (e.g., body size, ecological niche) associated with larger spikes?

- **Geography‑dependent evolution:**  
  Do lineages in region A exhibit larger spikes than those in region B?

---

## Dependencies

- **BDMM-Prime**  
  Vaughan & Stadler (2025). *Bayesian phylodynamic inference of multi-type population trajectories using genomic data.*  
  *Molecular Biology and Evolution* 42: msaf130  
  https://doi.org/10.1093/molbev/msaf130

- **BEAST 2.7**

---

## Citation

If you use this model, please cite:

- Manuscript in preparation (TBC)  
- Original Gamma Spike model:  
  Douglas, J., Bouckaert, R., Harris, S., Carter, C., & Wills, P. (2025).  
  *Evolution is coupled with branching across many granularities of life.*  
  Proc. R. Soc. B. 292:20250182  
  http://doi.org/10.1098/rspb.2025.0182
