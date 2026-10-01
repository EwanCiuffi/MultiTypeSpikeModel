# Multi-type Spike Model

**MultiTypeSpikeModel** is a BEAST 2 package for modelling **punctuated evolution**: bursts of change
("spikes") at branching events, alongside gradual change along branches under a relaxed clock. It
extends the Gamma Spike model of Douglas et al. (2025) so that:

- the size of spikes can depend on a lineage's **type** (e.g. a trait or a region), under the
  multi-type birth–death–migration tree prior of **BDMM-Prime**;
- birth, death and sampling rates can vary through time (**skyline** parameters);

The number and timing of hidden branching events (lineages that left no sampled descendants) are
integrated over rather than sampled.

---

## Model

The evolutionary distance along branch $e$ is

$$
d^e = r^e \mu_c \tau^e + \gamma^e \sum_{i=1}^{M} \mathbb{I}_i \, S^\mu_i \, \tilde{s}^e_i ,
$$

where $\mu_c$ is the mean clock rate, $r^e$ the relative branch rate, $\tau^e$ the branch duration,
$\gamma^e = 0$ for sampled-ancestor branches and 1 otherwise, $\mathbb{I}_i$ the spike indicator,
$S^\mu_i$ the spike mean and $\tilde{s}^e_i$ the latent spike of type $i$ ($M$ types).

Each branching event contributes a latent spike drawn from $\mathrm{Gamma}(S^\alpha_i, 1/S^\alpha_i)$.
With $N^e_i$ branching events of type $i$ on the branch (the observed event at its start, if of type
$i$, plus hidden events),

$$
\tilde{s}^e_i \sim \mathrm{Gamma}\left(N^e_i S^\alpha_i,\ \tfrac{1}{S^\alpha_i}\right).
$$

The number of hidden events is Poisson with a mean given by the tree prior, and the type of the
observed event is weighted by its probability under the tree prior; both are summed over exactly.

### Key parameters

| Parameter | Meaning |
|---|---|
| $S^\mu$ or $S^\mu_i$ | Spike mean: expected change per branching event. Shared or type-specific. |
| $S^\alpha$ or $S^\alpha_i$ | Spike shape: larger values give more uniform spikes. Shared or type-specific. |
| $\mathbb{I}$ or $\mathbb{I}_i$ | Spike indicator. When estimated, its posterior gives the support for punctuated evolution against a relaxed clock. |
| $\mu_c$, $\sigma_r$ | Mean clock rate and standard deviation of the log-normal relaxed clock. |

Type-specific spike means test whether punctuated change differs between types, for example between
lineages with different traits or in different regions.

---

## Installation

The package requires **BEAST 2.7.5** or later and installs **BDMM-Prime** and **SA** as dependencies.

### Through the BEAST Package Manager

1. In BEAUti, open **File > Manage Packages**, click **Package repositories** and **Add URL**, and enter
   ```
   https://raw.githubusercontent.com/EwanCiuffi/MultiTypeSpikeModel/main/packages.xml
   ```
2. Select **MultiTypeSpikeModel** in the package list and click **Install/Upgrade**.
3. Restart BEAUti.

### From source

You need OpenJDK 17 or later, the JavaFX SDK and Apache Ant. The build expects `beast2`, `BeastFX` and
`BDMM-Prime` as sibling directories of this repository, with BDMM-Prime built first:

```bash
git clone https://github.com/CompEvol/beast2.git
git clone https://github.com/CompEvol/BeastFX.git
git clone https://github.com/tgvaughan/BDMM-Prime.git
git clone https://github.com/EwanCiuffi/MultiTypeSpikeModel.git

(cd BDMM-Prime && ant)
cd MultiTypeSpikeModel
JAVA_FX_HOME=/path/to/javafx-sdk/lib ant install
```

`ant install` builds and tests the package, installs it into the BEAST 2.7 package directory and
resets BEAUti's cached package list; restart BEAUti afterwards. `ant package` builds the release zip
in `build/dist/` without installing it.

---

## Using the model in BEAUti

1. On the **Priors** tab, choose **BDMM-Prime** as the tree prior and set up its types and
   parameterization.
2. On the **Clock Model** tab, choose **MultiTypeSpikeClock (needs BDMM-Prime)**. The two steps can be
   done in either order; the spike prior always uses the BDMM-Prime settings of the same partition.
3. Adjust the priors on the spike mean, spike shape and clock parameters on the **Priors** tab.

Notes:

- The spike mean and spike shape are shared across types by default. To make them type-specific, set
  their dimension to the number of types in the XML.
- Estimate the indicator to compare the spike model with a relaxed clock without spikes.
- `nonCentered` switches the relaxed clock to a non-centred parameterisation, which can mix better when
  the branch rates are weakly informed by the data.

### Logged quantities

| Logger | Output |
|---|---|
| `SpikeLogger` | Spike on each branch, scaled by the spike mean (tree log). |
| `HiddenEventsLogger` | Number of hidden branching events per branch, sampled given the spikes, or their expectation with `logExpectedValue="true"` (tree log). |
| `SaltativeProportionLogger` | Proportion of the total distance due to spikes, overall and by type. |
| `NodeTypeProbabilityLogger` | Type probabilities of each node (multi-type analyses). |



## Citation

If you use this package, please cite:

- **Multi-type Spike Model**
  Ciuffi, E., Bickel, B., Vaughan, T. G., & Stadler, T. (in preparation).
  *Modelling trait-dependent punctuated evolution: the role of terrain ruggedness in Indo-European
diversification.*

- **Gamma Spike Model**
  Douglas, J., Bouckaert, R., Harris, S. C., Carter, C. W., & Wills, P. R. (2025).
  *Evolution is coupled with branching across many granularities of life.*
  Proceedings of the Royal Society B 292: 20250182.
  https://doi.org/10.1098/rspb.2025.0182

- **BDMM-Prime**
  Vaughan, T. G., & Stadler, T. (2025).
  *Bayesian phylodynamic inference of multi-type population trajectories using genomic data.*
  Molecular Biology and Evolution 42: msaf130.
  https://doi.org/10.1093/molbev/msaf130

- **BEAST 2**
  Bouckaert, R., Vaughan, T. G., Barido-Sottani, J., Duchêne, S., Fourment, M., Gavryushkina, A., ... & Drummond, A. J. (2019).
  *BEAST 2.5: An advanced software platform for Bayesian evolutionary analysis.*
  PLoS Computational Biology 15(4): e1006650.
  https://doi.org/10.1371/journal.pcbi.1006650
