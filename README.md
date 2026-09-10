# laplace-survival

A [Laplace](https://github.com/mlatinov/laplace) library of parametric survival models for Stan — Gompertz, log-logistic, Gompertz–Makeham, and piecewise-exponential — with right-censored likelihoods, pointwise log-likelihoods for LOO, survival curves, and event-time simulators. Import it into any `.laplace` model and call it with namespaced calls (`survival::function_name(...)`).

Like all Laplace libraries, `survival` compiles down to plain, readable Stan functions. Nothing about how you use it hides what actually ends up in your `.stan` file.

## Models

Every model is described by two functions: the **log hazard** $\log h(t)$, the instantaneous risk at time $t$, and the **cumulative hazard** $H(t)$, the total risk accumulated up to $t$. Everything else in the library — likelihoods, survival curves, simulation — is built from those two, using $S(t) = e^{-H(t)}$ and the right-censored log likelihood

$$
\log L = \sum_{i=1}^{N} \Big[ d_i \log h(t_i) - H(t_i) \Big],
$$

where $d_i = 1$ if the event was observed and $d_i = 0$ if the observation was right-censored.

| Model | Hazard $h(t)$ | Cumulative hazard $H(t)$ | Parameters |
|---|---|---|---|
| Gompertz | $e^{\eta + \gamma t}$ | $e^{\eta}\frac{e^{\gamma t}-1}{\gamma}$ | `gamma` log-hazard slope, `eta` log hazard at $t=0$ |
| Log-logistic | $\frac{(\beta/\alpha)(t/\alpha)^{\beta-1}}{1+(t/\alpha)^{\beta}}$ | $\log\left(1+(t/\alpha)^{\beta}\right)$ | `alpha` scale (median), `beta` shape |
| Gompertz–Makeham | $a + e^{\eta + \gamma t}$ | $at + e^{\eta}\frac{e^{\gamma t}-1}{\gamma}$ | `a` background hazard, `gamma`, `eta` as Gompertz |
| Piecewise exponential | $\lambda_{j} e^{\eta}$ for $t \in (s_j, s_{j+1}]$ | $e^{\eta}\sum_j \lambda_j E_j$ | `log_lambda` baseline log hazard per interval, `eta` linear predictor |

In the Gompertz, Makeham, and PWE models, `eta` is where your covariates go (e.g. `eta0 + X * beta`), giving a proportional-hazards model. The Gompertz and Makeham cumulative hazards are computed in a form that stays numerically stable as $\gamma \to 0$, where both reduce to an exponential model.

## What's included

All functions live in `survival_hazard.laplacelib` and follow the naming pattern `srv_<model>_<quantity>`, with `<model>` one of `gompertz`, `log_logistics`, `makeham`, `pwe`.

| Function | Returns | Gompertz | Log-logistic | Makeham | PWE |
|---|---|:-:|:-:|:-:|:-:|
| `srv_<model>_log_hazard` | $\log h(t)$ | ✓ | ✓ | ✓ | ✓ |
| `srv_<model>_cumulative_hazard` | $H(t)$ | ✓ | ✓ | ✓ | ✓ |
| `srv_<model>_lpdf` | Summed right-censored log likelihood | ✓ | ✓ | ✓ | `srv_pwe_loglik_sum` |
| `srv_<model>_lccdf` | Summed log survival, $\sum \log S(t_i)$ | ✓ | ✓ | ✓ | — |
| `srv_<model>_loglik` | Per-observation log likelihood (for LOO / WAIC) | ✓ | ✓ | ✓ | ✓ |
| `srv_<model>_survival_curve` | $S(t)$ | ✓ | ✓ | ✓ | ✓ |
| `srv_<model>_rng` | One simulated event time | ✓ | ✓ | ✓ | — |

The piecewise-exponential model adds two data-preparation helpers, meant to be called once in `transformed data`:

| Function | Returns |
|---|---|
| `srv_pwe_interval(t, starts)` | Interval index of each observation: the largest $j$ with $s_j < t$ |
| `srv_pwe_exposure(t, starts)` | $N \times J$ matrix of time spent at risk in each interval |

Most functions are overloaded so shared parameters can be passed as `real` and per-observation parameters as `vector`. Every function carries `@brief`, `@param`, `@return`, `@example`, and `@math` documentation, so you can read it from the terminal without leaving your model:

```
laplace doc survival::srv_gompertz_lpdf
```

### Argument order

Argument order is not identical across models — in particular, where the event indicator `d` goes — so here it is in one place:

| Model | Hazards / survival curve | `lpdf` | `lccdf` | `loglik` | `rng` |
|---|---|---|---|---|---|
| Gompertz | `(t, gamma, eta)` | `(t \| gamma, eta, d)` | `(t \| gamma, eta)` | `(t, d, gamma, eta)` | `(gamma, eta)` |
| Log-logistic | `(t, alpha, beta)` | `(t \| alpha, beta, d)` | `(t \| alpha, beta)` | `(t, d, alpha, beta)` | `(alpha, beta)` |
| Makeham | `(t, a, gamma, eta)` | `(t \| a, gamma, eta, d)` | `(t \| a, gamma, eta)` | `(t, a, gamma, eta, d)` | `(log_c, gamma, eta)` |
| PWE | `(interval, log_lambda, eta)` / `(E, log_lambda, eta)` | `loglik_sum(interval, E, d, log_lambda, eta)` | — | `(interval, E, d, log_lambda, eta)` | — |

## Installation

`survival` is distributed as a git-hosted Laplace library — there's no published registry entry yet, so it's added by pointing `laplace` (or `cmdlaplacer`, if you're working from R) directly at the repository. The package lives in the repository's `laplace/` subdirectory, so pass it as the subdir.

### Via the `laplace` CLI

From inside a Laplace project (a directory with its own `laplace.toml`):

```
laplace add survival --git https://github.com/mlatinov/laplace-survival --tag 0.1.0 --subdir laplace
```

### Via R (`cmdlaplacer`)

```r
library(cmdlaplacer)

laplace_install_git(
  "survival",
  "https://github.com/mlatinov/laplace-survival",
  tag = "0.1.0",
  subdir = "laplace"
)
```

Either way, this pins the dependency in your project's `laplace.toml`/`laplace.lock` at tag `0.1.0`. Check the [releases](https://github.com/mlatinov/laplace-survival/tags) for newer tags as they become available.

## Usage

Import the library in a `library { }` block and call its functions with the `survival::` namespace prefix.

### Gompertz proportional-hazards model

A Gompertz model with covariates, right censoring, pointwise log likelihoods for `loo`, a baseline survival curve, and a simulated event time:

```stan
library {
    import survival
}

data {
  int<lower=1> N;                      // subjects
  int<lower=1> K;                      // covariates
  vector<lower=0>[N] t;                // follow-up time
  vector<lower=0, upper=1>[N] d;       // 1 = event, 0 = censored
  matrix[N, K] X;
  int<lower=1> G;
  vector<lower=0>[G] t_grid;           // times for the survival curve
}

parameters {
  real eta0;                           // baseline log hazard at t = 0
  vector[K] beta;                      // log hazard ratios
  real gamma;                          // change in log hazard per unit time
}

model {
  vector[N] eta = eta0 + X * beta;

  eta0 ~ normal(-3, 2);
  beta ~ normal(0, 1);
  gamma ~ normal(0, 0.5);

  target += survival::srv_gompertz_lpdf(t | gamma, eta, d);
}

generated quantities {
  vector[N] log_lik = survival::srv_gompertz_loglik(t, d, gamma, eta0 + X * beta);
  vector[G] S_baseline = survival::srv_gompertz_survival_curve(t_grid, gamma, rep_vector(eta0, G));
  real t_new = survival::srv_gompertz_rng(gamma, eta0);
}
```

Swapping in another parametric model is a matter of changing the function family and its parameters, e.g. `survival::srv_log_logistics_lpdf(t | alpha, beta, d)` or `survival::srv_makeham_lpdf(t | a, gamma, eta, d)`.

### Piecewise-exponential model

The interval lookup and exposure matrix depend only on data, so they are computed once in `transformed data`:

```stan
library {
    import survival
}

data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> J;                      // number of intervals
  vector<lower=0>[N] t;
  vector<lower=0, upper=1>[N] d;
  matrix[N, K] X;
  vector<lower=0>[J] starts;           // interval start points, starts[1] = 0
}

transformed data {
  array[N] int interval = survival::srv_pwe_interval(t, starts);
  matrix[N, J] E = survival::srv_pwe_exposure(t, starts);
}

parameters {
  vector[J] log_lambda;                // baseline log hazard per interval
  vector[K] beta;
}

model {
  log_lambda ~ normal(-3, 1);
  beta ~ normal(0, 1);

  target += survival::srv_pwe_loglik_sum(interval, E, d, log_lambda, X * beta);
}

generated quantities {
  vector[N] log_lik = survival::srv_pwe_loglik(interval, E, d, log_lambda, X * beta);
}
```

### From R

With `cmdlaplacer`, the `.laplace` file compiles straight to a `cmdstanr` model, and the generated `.stan` file stays on disk next to it:

```r
library(cmdlaplacer)

mod <- laplace_model("gompertz.laplace")
fit <- mod$sample(data = stan_data)

fit$loo()   # uses the log_lik vector from generated quantities
```

## Things to know

- **Use `target +=`, not `~`.** Laplace rewrites `survival::func(` calls, so write `target += survival::srv_gompertz_lpdf(t | ...)`. The `t ~ survival::srv_gompertz(...)` form won't resolve.
- **The `lpdf` functions already handle censoring** through `d`, so pass every observation — events and censored — in one call. The `lccdf` functions are for the alternative pattern where censored rows are passed separately; never use both for the same rows, or censored observations are counted twice.
- **`d` is a `vector`**, not an integer array. If your indicator is `array[N] int`, convert it once in `transformed data` with `to_vector(...)`.
- **Right censoring only.** Left censoring, interval censoring, and delayed entry (left truncation) are not supported yet.
- **Log-logistic times must be strictly positive**, since the hazard involves $\log t$.
- **Gompertz with $\gamma < 0$ implies a cure fraction:** some subjects never experience the event, and `srv_gompertz_rng` returns `positive_infinity()` for them. Account for this when summarising simulated times.
- **`srv_makeham_rng` takes the background hazard on the log scale** (`log_c`), unlike the other Makeham functions, which take `a` on the natural scale. Pass `log(a)`.
- **PWE intervals** should start at 0 (`starts[1] = 0`). A time exactly on a boundary belongs to the earlier interval, and the last interval is open-ended.
- **Survival curves take one `eta` per time point.** For a single covariate profile over a time grid, use `rep_vector(eta0, G)`.

## License

See [LICENSE](./LICENSE).