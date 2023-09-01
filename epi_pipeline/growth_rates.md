Let $x_d$ be the number of sequences of a variant in a location on date $d$ and $b_d$ be the number of non-variant sequences. Let $n_d$ be the estimated total number of new coronavirus cases in that location on that date; then we can estimate the number of new cases of the query variant $x'_d$ and the associated uncertainty by treating each unsequenced case as a Bernoulli trial. This results in a binomial distribution, which can be closely approximated by a normal distribution $\mathcal{N}(\frac{n_d x_d}{x_d + b_d}, \frac{\sqrt{n_d x_d b_d}}{x_d + b_d})$. We can then use this model with uncertainty propagation to estimate the uncertainty in $\text{ln}(x_d)$ as $(b_d / x_d^3)^{0.25}$. We then select weights for time-wise smoothing based on 

We can then use this model to estimate the uncertainty in $x_d$ (or $b_d$) as $\sqrt{\frac{x_d b_d}{n_d}}$. 

$x_d = x'_d (x_d+b_d) / n_d$
sqrt
$\mathcal{N}(x_d, \frac{x_d \sqrt{b_d}}{x_d + b_d})$


$\mathcal{N}(x_d, \frac{x_d \sqrt{b_d}}{x_d + b_d})$




$\frac{n_d x_d}{x_d + b_d}$ is , and we can estimate $\Delta \frac{n_d x_d}{x_d + b_d}$ by . 
