# The Latent Class Model

## Latent Class Analysis

Latent Class Analysis assumes that people belong to different groups
(i.e., classes) that are due to nonobservable characteristics. In a
latent class model, we aim to estimate two kinds of probabilities:

1.  The probability that a person belongs to a particular class.
2.  The conditional probabilities of item responses if a person belongs
    to a given class.

## The likelihood

Suppose that a sample of people respond to $J$ items and $\mathbf{y}$ is
a vector that contains the scores to each item $j$. Let $K$ denote the
number of latent classes and let $x_{k}$ be class $k$. Then, the
likelihood of the response pattern $\mathbf{y}$, if it was observed $n$
times in the sample and every person responds independently of each
other, can be written as
$$\ell\; = \; P(\mathbf{y})^{n}\; = \;\left( \sum\limits_{k = 1}^{K}P(x_{k})\, P(\mathbf{y} \mid x_{k}) \right)^{n}.$$

Assuming **local independence** (responses within a person are
independent conditional on class), we have
$$P(\mathbf{y} \mid x_{k})\; = \;\prod\limits_{j = 1}^{J}P(y_{j} \mid x_{k}),$$
where $y_{j}$ denotes the score on item $j$. Hence,
$$\ell\; = \;\left( \sum\limits_{k = 1}^{K}P(x_{k})\,\prod\limits_{j = 1}^{J}P(y_{j} \mid x_{k}) \right)^{\! n},$$
and the log-likelihood becomes
$$\ell\ell\; = \; n{\log}\left( \sum\limits_{k = 1}^{K}P(x_{k})\,\prod\limits_{j = 1}^{J}P(y_{j} \mid x_{k}) \right).$$
The term inside the parenthesis is the probability of a single pattern,
$P(\mathbf{y})$. Assuming independence between people with different
response patterns, the log-likelihood of the whole sample is the sum of
log-likelihoods of each response pattern.

To simplify the computation of the logarithm likelihood and related
derivatives, let $l_{y_{m} \mid x_{g}} = {\log}P(y_{m} \mid x_{g})$, so
that

$$\ell\ell\; = \; n{\log}\left( \sum\limits_{k = 1}^{K}P(x_{k})\,{\exp}\left( \sum\limits_{j = 1}^{J}l_{y_{m} \mid x_{g}} \right) \right).$$

## First-order derivatives

For a fixed pattern $\mathbf{y}$, define $$\begin{aligned}
{P(\mathbf{y})\;} & {= \;\sum\limits_{k = 1}^{K}P(x_{k})\,\prod\limits_{j = 1}^{J}P(y_{j} \mid x_{k})} \\
 & {= \;\sum\limits_{k = 1}^{K}P(x_{k})\,{\exp}\left( \sum\limits_{j = 1}^{J}l_{y_{m} \mid x_{g}} \right).}
\end{aligned}$$ Then
$$\frac{\partial\ell\ell}{\partial P(x_{g})}\; = \; n\,\frac{1}{P(\mathbf{y})}\,\prod\limits_{j = 1}^{J}P(y_{j} \mid x_{g}).$$
For a specific item $m$ and class $g$,
$$\frac{\partial\ell\ell}{\partial l_{y_{m} \mid x_{g}}}\; = \; n\,\frac{1}{P(\mathbf{y})}\, P(x_{g})\,\prod\limits_{j = 1}^{J}P(y_{j} \mid x_{g}).$$

Notice that this last expression is just the posterior,
$P(x_{g} \mid y_{m})$, weighted by $n$.

## Directional derivatives of first-order derivatives

$$d\left\lbrack \frac{\partial\ell\ell}{\partial P(x_{g})} \right\rbrack\; = \; n\,\frac{P(\mathbf{y})\; d\left\lbrack \prod\limits_{j = 1}^{J}P(y_{j} \mid x_{g}) \right\rbrack - d\left\lbrack P(\mathbf{y}) \right\rbrack\;\prod\limits_{j = 1}^{J}P(y_{j} \mid x_{g})}{P(\mathbf{y})^{2}},$$
where

$$d\left\lbrack \prod\limits_{j = 1}^{J}P(y_{j} \mid x_{g}) \right\rbrack = \frac{\prod\limits_{j = 1}^{J}P(y_{j} \mid x_{g})}{P(y_{j} \mid x_{g})\, d\left\lbrack P(y_{j} \mid x_{g}) \right\rbrack},$$
and

$$d\left\lbrack P(\mathbf{y}) \right\rbrack = \sum\limits_{k = 1}^{K}d\left\lbrack P(x_{k})\rbrack\,\prod\limits_{j = 1}^{J}P(y_{j} \mid x_{k}) + \sum\limits_{k = 1}^{K}P(x_{k})\, d\left\lbrack \prod\limits_{j = 1}^{J}P(y_{j} \mid x_{k})\rbrack. \right. \right.$$

## Second-order derivatives

For classes $g,h$, $$\begin{aligned}
{\frac{\partial^{2}\ell\ell}{\partial P(x_{g})\,\partial P(x_{h})}\;} & {= \; - \, n\,\frac{1}{P(\mathbf{y})^{2}}\,\left( \prod\limits_{j = 1}^{J}P(y_{j} \mid x_{g}) \right)\!\left( \prod\limits_{j = 1}^{J}P(y_{j} \mid x_{h}) \right)} \\
 & {= \; - n\frac{P(x_{g} \mid y)P(x_{h} \mid y)}{P(x_{g})P(x_{h})}.}
\end{aligned}$$

For items $m,n$ and classes $g,h$,
$$\frac{\partial^{2}\ell\ell}{\partial l_{y_{m} \mid x_{g}}\,\partial l_{y_{n} \mid x_{h}}} = \begin{cases}
{\, n\, P(x_{g} \mid y)\!\ \left( 1 - P(x_{g} \mid y) \right),} & {{\text{if}\mspace{6mu}}\!\ g = h,} \\
{-\, n\, P(x_{g} \mid y)\!\ P(x_{h} \mid y),} & \text{otherwise.}
\end{cases}$$

For the mixed second derivative,
$$\frac{\partial^{2}\ell\ell}{\partial P(x_{h})\,\partial l_{y_{m} \mid x_{g}}} = \begin{cases}
{\frac{1}{P(x_{g})}\, n\, P(x_{g} \mid y)\left( 1 - P(x_{g} \mid y) \right),} & {{\text{if}\mspace{6mu}}g = h,} \\
{-\,\frac{1}{P(x_{h})}\, n\, P(x_{g} \mid y)P(x_{h} \mid y),} & \text{otherwise.}
\end{cases}$$

Collecting these terms gives the Hessian in block form:
$${Hess}(\ell\ell)\; = \;\begin{bmatrix}
\frac{\partial^{2}\ell\ell}{\partial P(x_{g})\,\partial P(x_{h})} & \frac{\partial^{2}\ell\ell}{\partial P(x_{h})\,\partial l_{y_{m} \mid x_{g}}} \\
\frac{\partial^{2}\ell\ell}{\partial l_{y_{m} \mid x_{g}}\,\partial P(x_{h})} & \frac{\partial^{2}\ell\ell}{\partial l_{y_{m} \mid x_{g}}\,\partial l_{y_{m} \mid x_{g}}}
\end{bmatrix}\!.$$

## Models for the conditional likelihoods

The conditional probabilities need to be parameterized with a
likelihood. We consider a **multinomial** likelihood for categorical
items and a **Gaussian** likelihood for continuous items.

### Multinomial

For categorical items, let $\pi_{m_{k} \mid g}$ be the probability of
scoring category $k$ on item $m$ if a subject belongs to class $g$. Then
$$P(y_{m} \mid x_{g})\; = \;\pi_{m_{k} \mid g},$$ where $k$ is such that
$y_{m} = k$. With this parameterization,
$$\frac{\partial l_{y_{m} \mid x_{g}}}{\partial\pi_{n_{k} \mid h}} = \begin{cases}
{\frac{1}{\pi_{n_{k} \mid h}},} & {{\text{if}\mspace{6mu}}y_{m} = k,\ m = n,\ g = h,} \\
{0,} & \text{otherwise.}
\end{cases}$$ and
$$\frac{\partial^{2}l_{y_{j} \mid x_{i}}}{\partial\pi_{m_{k} \mid g}\partial\pi_{n_{l} \mid h}} = \begin{cases}
{-\frac{1}{\pi_{m_{j} \mid i}^{2}},} & {{\text{if}\mspace{6mu}}y_{j} = k = l,\ y = m = n,\ i = g = h,} \\
{0,} & \text{otherwise.}
\end{cases}$$

Consequently, the Hessian for each conditional parameter has the
following block form:
$${Hess}(l_{y_{m} \mid x_{g}})\; = \;\mathbf{e}\;\frac{\partial^{2}l_{y_{j} \mid x_{i}}}{\partial\pi_{m_{k} \mid g}\partial\pi_{n_{l} \mid h}}\;\mathbf{e}^{\top},$$
where $\mathbf{e}$ is a vector of zeroes with a 1 in the position
corresponding to the parameter $\pi_{y_{j} \mid i}$.

Notice that each conditional parameter $l_{y_{m} \mid x_{g}}$ has a
Hessian matrix.

### Gaussian

For continuous items, let $\varphi$ denote the normal density. Let
$\mu_{m \mid g}$ and $\sigma_{m \mid g}$ be the mean and standard
deviation for item $m$ in class $g$. Then
$$P(y_{m} \mid x_{g})\; = \;\varphi\!\left( y_{m};\,\mu_{m \mid g},\,\sigma_{m \mid g} \right).$$

First derivatives:
$$\frac{\partial l_{y_{m} \mid x_{g}}}{\partial\mu_{n \mid h}} = \begin{cases}
{\frac{y_{m} - \mu_{m \mid g}}{\sigma_{m \mid g}^{2}},} & {{\text{if}\mspace{6mu}}m = n,\ g = h,} \\
{0,} & \text{otherwise,}
\end{cases}$$

$$\frac{\partial l_{y_{m} \mid x_{g}}}{\partial\sigma_{n \mid h}} = \begin{cases}
{\frac{(y_{m} - \mu_{m \mid g})^{2} - \sigma_{m \mid g}^{2}}{\sigma_{m \mid g}^{3}},} & {{\text{if}\mspace{6mu}}m = n,\ g = h,} \\
{0,} & \text{otherwise.}
\end{cases}$$

Second-order derivatives:
$$\frac{\partial^{2}l_{y_{m} \mid x_{g}}}{\partial\mu_{m \mid g}\,\partial\mu_{n \mid h}} = \begin{cases}
{-\frac{1}{\sigma_{m \mid g}^{2}},} & {{\text{if}\mspace{6mu}}m = n,\ g = h,} \\
{0,} & \text{otherwise,}
\end{cases}$$

$$\frac{\partial^{2}l_{y_{m} \mid x_{g}}}{\partial\sigma_{m \mid g}\,\partial\sigma_{n \mid h}} = \begin{cases}
{\frac{1}{\sigma_{m \mid g}^{2}} - \frac{3(y_{m} - \mu_{m \mid g})^{2}}{\sigma_{m \mid g}^{4}},} & {{\text{if}\mspace{6mu}}m = n,\ g = h,} \\
{0,} & \text{otherwise,}
\end{cases}$$

$$\frac{\partial^{2}l_{y_{m} \mid x_{g}}}{\partial\mu_{m \mid g}\,\partial\sigma_{n \mid h}} = \begin{cases}
{-\frac{2(y_{m} - \mu_{m \mid g})}{\sigma_{m \mid g}^{3}},} & {{\text{if}\mspace{6mu}}m = n,\ g = h,} \\
{0,} & \text{otherwise.}
\end{cases}$$

Consequently, the Hessian for each conditional parameter has the
following block form:
$${Hess}(l_{y_{m} \mid x_{g}})\; = \;\begin{bmatrix}
\frac{\partial^{2}l_{y_{m} \mid x_{g}}}{\partial\mu_{m \mid g}\,\partial\mu_{n \mid h}} & \frac{\partial^{2}l_{y_{m} \mid x_{g}}}{\partial\mu_{m \mid g}\,\partial\sigma_{n \mid h}} \\
\frac{\partial^{2}l_{y_{m} \mid x_{g}}}{\partial\sigma_{n \mid h}\,\partial\mu_{m \mid g}} & \frac{\partial^{2}l_{y_{m} \mid x_{g}}}{\partial\sigma_{m \mid g}\,\partial\sigma_{n \mid h}}
\end{bmatrix}.$$

Notice that each conditional parameter $l_{y_{m} \mid x_{g}}$ has a
Hessian matrix.

## Model for the latent class probabilities

Probabilities of class membership are parameterized with the softmax
transformation:

$$P(x_{g}) = \frac{{\exp}(\theta_{g})}{\sum\limits_{j}{\exp}(\theta_{j})},$$
where $\theta_{g}$ is the log-scale parameter associated with class $g$.

The jacobian of this transformation is given by

$$J = {diag}(P) - PP^{\top}.$$ Finally, the Hessian for each probability
is

$${Hess}\left( P(x_{g}) \right)\; = \; P(x_{g})\left( (e_{g} - P)(e_{g} - P)^{\top} - J \right),$$
where $e_{g}$ is a vector of zeroes with a $1$ in position $g$.

## Model for the conditional probabilities of the multinomial model

Probabilities of conditional responses are parameterized with the
softmax transformation:

$$\pi_{m_{k} \mid g} = \frac{{\exp}(\eta_{m_{k} \mid g})}{\sum\limits_{j}{\exp}(\eta_{m_{j} \mid g})},$$
where $\eta_{m_{k} \mid g}$ is the log-scale parameter associated with
response $k$ to item $m$ in class $g$.

The jacobian of this transformation is given by

$$J_{m \mid g} = {diag}(\pi_{m \mid g}) - \pi_{m \mid g}(\pi_{m \mid g})^{\top}.$$
Finally, the Hessian for each probability is

$${Hess}\left( \pi_{m_{k} \mid g} \right)\; = \;\pi_{m \mid g}\left( (e_{k} - \pi_{m \mid g})(e_{k} - \pi_{m \mid g})^{\top} - J \right),$$
where $e_{k}$ is a vector of zeroes with a $1$ in position $k$.

## Constant priors

### For latent class probabilities

### For conditional likelihoods

#### Multinomial

For the conditional probabilities modeled with a multinomial likelihood,
we add the following term to the log-likelihood for each class $g$:

$$\lambda_{2} = \sum\limits_{m}\sum\limits_{k}{\widehat{\pi}}_{m_{k}}\sum\limits_{g}\frac{\alpha}{K}{\log}(\pi_{m_{k} \mid g}),$$

Where ${\widehat{\pi}}_{m_{k}}$ is the proportion of times category $k$
was selected in item $m$.

The first-order derivatives are

$$\frac{\partial\lambda_{2}}{\partial\pi_{m_{k} \mid g}} = \sum\limits_{m}\sum\limits_{k}{\widehat{\pi}}_{m_{k}}\sum\limits_{g}\frac{\alpha}{K}\frac{1}{\pi_{m_{k} \mid g}}.$$

The second-order derivatives are

$$\frac{\partial^{2}\lambda_{2}}{\partial\pi_{m_{k} \mid g}\;\partial\pi_{m_{k} \mid g}} = -\sum\limits_{m}\sum\limits_{k}{\widehat{\pi}}_{m_{k}}\sum\limits_{g}\frac{\alpha}{K}\frac{1}{\pi_{m_{k} \mid g}^{2}}.$$

#### Gaussian

For the conditional probabilities modeled with a gaussian likelihood, we
add the following term to the log-likelihood for each class $g$:
$$\begin{aligned}
\lambda_{3} & {= \sum\limits_{g}^{K}\left( -0.5\frac{\alpha}{K}\log\left( \prod\limits_{j}\sigma_{j \mid g}^{2} \right) - 0.5\frac{\alpha}{K}\sum\limits_{j}\frac{{\widehat{\sigma}}_{j}^{2}}{\sigma_{j \mid g}^{2}} \right)} \\
 & {= \sum\limits_{g}^{K}\left( -\frac{\alpha}{K}\sum\limits_{j}s_{j \mid g} - 0.5\frac{\alpha}{K}\sum\limits_{j}\frac{{\widehat{\sigma}}_{j}^{2}}{\sigma_{j \mid g}^{2}} \right),}
\end{aligned}$$

where ${\widehat{\sigma}}_{j}^{2}$ is the variance of item $j$.

The first-order derivatives are $$\begin{aligned}
{\frac{\partial\lambda_{3}}{\partial\sigma_{m \mid g}}\;} & {= \;\sum\limits_{g}^{K} - \frac{\alpha}{K\sigma_{m \mid g}} + \frac{\alpha}{K\sigma_{m \mid g}^{3}}{\widehat{\sigma}}_{m \mid g}^{2}} \\
 & {= \sum\limits_{g}^{K}\frac{\alpha}{K}\left( \frac{{\widehat{\sigma}}_{m \mid g}^{2}}{\sigma_{m \mid g}^{3}} - \frac{1}{\sigma_{m \mid g}} \right).}
\end{aligned}$$

The second-order derivatives are
$$\frac{\partial^{2}\lambda_{3}}{\partial\sigma_{m \mid g}\partial\sigma_{m \mid g}}\; = \;\sum\limits_{g}^{K}\frac{\alpha}{K}\left( \frac{1}{\sigma_{m \mid g}} - 3\frac{{\widehat{\sigma}}_{m \mid g}^{2}}{\sigma_{m \mid g}^{4}} \right).$$
