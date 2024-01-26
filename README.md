# EccentricAmplitudes
Eccentric amplitudes for an inspiraling black hole binary.
Follows Arredondo _et al._ (in prep.).

## Requirements and usage
This code was prepared with Mathematica 13.

Make sure this directory's parent is in your notebook's path with ```AppendTo[$Path, 
  "/path/to/parent/directory"]```, or place this directory
  within your Mathematica's contents.
One can then load the package with ```<< EccentricAmplitudes` ```.

## Package contents
This package is designed to easily recover the harmonic decomposition of the reduced GW modes,
$$
\hat{H}^{lm} = \sum_{j=-j_{\rm min}}^{j_{\rm max}} N^{lm}_j \exp(-ij\ell).
$$
The Fourier amplitudes $N^{lm}_j$ are obtained in table with the function `AmplitudesTable` for a given eccentricity and can be returned symbolically or evaluated numerically.
Along with the table, a list of the harmonics $j$ is returned with it for easy identification of each term.
Use `?AmplitudesTable` to see its usage.
We recommend numerical evaluation if speed is of the essence.

The full series itself (as written above) can be evaluated with `recoverSeries`.

### The nested limits
The series from above are enough to reconstruct the GW modes, but one can also dive into the details of the Fourier series.
The full list of series limits within each $j$ can be found with `recoverConvList`.
The list for $m=0$ is simple, but for $m=2$ is intricate due to the nested limits involved.
The structure for the latter is as follows, given an eccentricity:

{{ $j$ limits }, { $s$ limits for each $j$ }, { $k$ limits for $\mathcal{P}^{mW}_s$ and $\mathcal{K}^{lm}_{s-j}$ at each $(s,j)$ pair }}

For example, for $e=0.01$,

```mathematica
 In:= recoverConvList[1/100,2,2]
Out:= {{-1,1},{{-2,2},{-1,2},{-1,3}},{{{2,2},{3,1},{4,1},{4,1},{5,1}},{{3,2},{4,1},{4,1},{5,1}},{{3,2},{4,2},{4,1},{5,1},{6,1}}}}
```
The three lists are returned in the output, of which we look at certain entries from each:
- Limits in $j$:
    - The top-level series ranges $j \in [-1,1]$
- Limits in $s$:
    - For $j=-1$, the series expansion $N^{22}_{-1} = \sum_s \mathcal{P}^{2W}_s \mathcal{K}^{22}_{s-j}$ ranges $s \in [-2,2]$
    - For $j=0$ it ranges $s \in [-1,2]$
- Limits in $k$:
    - For $j=-1,s=-2$: $\mathcal{P}^{2W}_{-2} = \sum_{k=0}^{k_{\rm max}} \mathcal{P}^{2W}_{-1,k}$ is truncated at $k_{\rm max} = 2$ while $\mathcal{K}^{22}_{-1} = \sum_{k=0}^{k_{\rm max}} \mathcal{K}^{22}_{-1,k}$ is truncated at $k_{\rm max}=2$ as well
    - For $j=-1,s=-1$: $\mathcal{P}^{2W}_{-1} = \sum_{k=0}^{k_{\rm max}} \mathcal{P}^{2W}_{-1,k}$ is truncated at $k_{\rm max} = 3$ while $\mathcal{K}^{22}_{0} = \sum_{k=0}^{k_{\rm max}} \mathcal{K}^{22}_{0,k}$ is truncated at $k_{\rm max}=1$
    - For $j=0,s=-1$: $\mathcal{P}^{2W}_{-1}$ is truncated at $k_{\rm max} = 3$ while $\mathcal{K}^{22}_{-1}$ is truncated at $k_{\rm max}=2$. Note that these specific terms appeared in the entries that we discussed above

