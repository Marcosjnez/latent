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
$$\ell\; = \; P(\mathbf{y})^{n}\; = \;(\sum\limits_{k = 1}^{K}P\left( x_{k} \right)\, P\left( \mathbf{y} \mid x_{k} \right))^{n}.$$

Assuming **local independence** (responses within a person are
independent conditional on class), we have
$$P\left( \mathbf{y} \mid x_{k} \right)\; = \;\prod\limits_{j = 1}^{J}P\left( y_{j} \mid x_{k} \right),$$
where $y_{j}$ denotes the score on item $j$. Hence,
$$\ell\; = \;(\sum\limits_{k = 1}^{K}P\left( x_{k} \right)\,\prod\limits_{j = 1}^{J}P\left( y_{j} \mid x_{k} \right))^{\! n},$$
and the log-likelihood becomes
$$\ell\ell\; = \; n\log(\sum\limits_{k = 1}^{K}P\left( x_{k} \right)\,\prod\limits_{j = 1}^{J}P\left( y_{j} \mid x_{k} \right)).$$
The term inside the parenthesis is the probability of a single pattern,
$P(\mathbf{y})$. Assuming independence between people with different
response patterns, the log-likelihood of the whole sample is the sum of
log-likelihoods of each response pattern.

To simplify the computation of the logarithm likelihood and related
derivatives, let
$l_{y_{m} \mid x_{g}} = \log P\left( y_{m} \mid x_{g} \right)$, so that

$$\ell\ell\; = \; n\log(\sum\limits_{k = 1}^{K}P\left( x_{k} \right)\,\exp(\sum\limits_{j = 1}^{J}l_{y_{m} \mid x_{g}})).$$

## First-order derivatives

For a fixed pattern $\mathbf{y}$, define $$\begin{aligned}
{P(\mathbf{y})\;} & {= \;\sum\limits_{k = 1}^{K}P\left( x_{k} \right)\,\prod\limits_{j = 1}^{J}P\left( y_{j} \mid x_{k} \right)} \\
 & {= \;\sum\limits_{k = 1}^{K}P\left( x_{k} \right)\,\exp(\sum\limits_{j = 1}^{J}l_{y_{m} \mid x_{g}}).}
\end{aligned}$$ Then
$$\frac{\partial\ell\ell}{\partial P\left( x_{g} \right)}\; = \; n\,\frac{1}{P(\mathbf{y})}\,\prod\limits_{j = 1}^{J}P\left( y_{j} \mid x_{g} \right).$$
For a specific item $m$ and class $g$,
$$\frac{\partial\ell\ell}{\partial l_{y_{m} \mid x_{g}}}\; = \; n\,\frac{1}{P(\mathbf{y})}\, P\left( x_{g} \right)\,\prod\limits_{j = 1}^{J}P\left( y_{j} \mid x_{g} \right).$$

Notice that this last expression is just the posterior,
$P\left( x_{g} \mid y_{m} \right)$, weighted by $n$.

## Directional derivatives of first-order derivatives

$$d\lbrack\frac{\partial\ell\ell}{\partial P\left( x_{g} \right)}\rbrack\; = \; n\,\frac{P(\mathbf{y})\; d\lbrack\prod\limits_{j = 1}^{J}P\left( y_{j} \mid x_{g} \right)\rbrack - d\lbrack P(\mathbf{y})\rbrack\;\prod\limits_{j = 1}^{J}P\left( y_{j} \mid x_{g} \right)}{P(\mathbf{y})^{2}},$$
where

$$d\lbrack\prod\limits_{j = 1}^{J}P\left( y_{j} \mid x_{g} \right)\rbrack = \frac{\prod\limits_{j = 1}^{J}P\left( y_{j} \mid x_{g} \right)}{P\left( y_{j} \mid x_{g} \right)\, d\lbrack P\left( y_{j} \mid x_{g} \right)\rbrack},$$
and

$$d\lbrack P(\mathbf{y})\rbrack = \sum\limits_{k = 1}^{K}d\lbrack P\left( x_{k} \right)\rbrack\,\prod\limits_{j = 1}^{J}P\left( y_{j} \mid x_{k} \right) + \sum\limits_{k = 1}^{K}P\left( x_{k} \right)\, d\lbrack\prod\limits_{j = 1}^{J}P\left( y_{j} \mid x_{k} \right)\rbrack.$$

## Second-order derivatives

For classes $g,h$, $$\begin{aligned}
{\frac{\partial^{2}\ell\ell}{\partial P\left( x_{g} \right)\,\partial P\left( x_{h} \right)}\;} & {= \; - \, n\,\frac{1}{P(\mathbf{y})^{2}}\,(\prod\limits_{j = 1}^{J}P\left( y_{j} \mid x_{g} \right))\!(\prod\limits_{j = 1}^{J}P\left( y_{j} \mid x_{h} \right))} \\
 & {= \; - n\frac{P\left( x_{g} \mid y \right)P\left( x_{h} \mid y \right)}{P\left( x_{g} \right)P\left( x_{h} \right)}.}
\end{aligned}$$

For items $m,n$ and classes $g,h$,
$$\frac{\partial^{2}\ell\ell}{\partial l_{y_{m} \mid x_{g}}\,\partial l_{y_{n} \mid x_{h}}} = \begin{cases}
{\, n\, P\left( x_{g} \mid y \right)\!\ (1 - P\left( x_{g} \mid y \right)),} & {{\text{if}\mspace{6mu}}\!\ g = h,} \\
{-\, n\, P\left( x_{g} \mid y \right)\!\ P\left( x_{h} \mid y \right),} & \text{otherwise.}
\end{cases}$$

For the mixed second derivative,
$$\frac{\partial^{2}\ell\ell}{\partial P\left( x_{h} \right)\,\partial l_{y_{m} \mid x_{g}}} = \begin{cases}
{\frac{1}{P\left( x_{g} \right)}\, n\, P\left( x_{g} \mid y \right)(1 - P\left( x_{g} \mid y \right)),} & {{\text{if}\mspace{6mu}}g = h,} \\
{-\,\frac{1}{P\left( x_{h} \right)}\, n\, P\left( x_{g} \mid y \right)P\left( x_{h} \mid y \right),} & \text{otherwise.}
\end{cases}$$

Collecting these terms gives the Hessian in block form:
$${Hess}\left( \ell\ell \right)\; = \;\begin{bmatrix}
\frac{\partial^{2}\ell\ell}{\partial P\left( x_{g} \right)\,\partial P\left( x_{h} \right)} & \frac{\partial^{2}\ell\ell}{\partial P\left( x_{h} \right)\,\partial l_{y_{m} \mid x_{g}}} \\
\frac{\partial^{2}\ell\ell}{\partial l_{y_{m} \mid x_{g}}\,\partial P\left( x_{h} \right)} & \frac{\partial^{2}\ell\ell}{\partial l_{y_{m} \mid x_{g}}\,\partial l_{y_{m} \mid x_{g}}}
\end{bmatrix}\!.$$

## Models for the conditional likelihoods

The conditional probabilities need to be parameterized with a
likelihood. We consider a **multinomial** likelihood for categorical
items and a **Gaussian** likelihood for continuous items.

### Multinomial

For categorical items, let $\pi_{m_{k} \mid g}$ be the probability of
scoring category $k$ on item $m$ if a subject belongs to class $g$. Then
$$P\left( y_{m} \mid x_{g} \right)\; = \;\pi_{m_{k} \mid g},$$ where $k$
is such that $y_{m} = k$. With this parameterization,
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
$${Hess}\left( l_{y_{m} \mid x_{g}} \right)\; = \;\mathbf{e}\;\frac{\partial^{2}l_{y_{j} \mid x_{i}}}{\partial\pi_{m_{k} \mid g}\partial\pi_{n_{l} \mid h}}\;\mathbf{e}^{\top},$$
where $\mathbf{e}$ is a vector of zeroes with a 1 in the position
corresponding to the parameter $\pi_{y_{j} \mid i}$.

Notice that each conditional parameter $l_{y_{m} \mid x_{g}}$ has a
Hessian matrix.

### Gaussian

For continuous items, let $\varphi$ denote the normal density. Let
$\mu_{m \mid g}$ and $\sigma_{m \mid g}$ be the mean and standard
deviation for item $m$ in class $g$. Then
$$P\left( y_{m} \mid x_{g} \right)\; = \;\varphi\!(y_{m};\,\mu_{m \mid g},\,\sigma_{m \mid g}).$$

First derivatives:
$$\frac{\partial l_{y_{m} \mid x_{g}}}{\partial\mu_{n \mid h}} = \begin{cases}
{\frac{y_{m} - \mu_{m \mid g}}{\sigma_{m \mid g}^{2}},} & {{\text{if}\mspace{6mu}}m = n,\ g = h,} \\
{0,} & \text{otherwise,}
\end{cases}$$

$$\frac{\partial l_{y_{m} \mid x_{g}}}{\partial\sigma_{n \mid h}} = \begin{cases}
{\frac{\left( y_{m} - \mu_{m \mid g} \right)^{2} - \sigma_{m \mid g}^{2}}{\sigma_{m \mid g}^{3}},} & {{\text{if}\mspace{6mu}}m = n,\ g = h,} \\
{0,} & \text{otherwise.}
\end{cases}$$

Second-order derivatives:
$$\frac{\partial^{2}l_{y_{m} \mid x_{g}}}{\partial\mu_{m \mid g}\,\partial\mu_{n \mid h}} = \begin{cases}
{-\frac{1}{\sigma_{m \mid g}^{2}},} & {{\text{if}\mspace{6mu}}m = n,\ g = h,} \\
{0,} & \text{otherwise,}
\end{cases}$$

$$\frac{\partial^{2}l_{y_{m} \mid x_{g}}}{\partial\sigma_{m \mid g}\,\partial\sigma_{n \mid h}} = \begin{cases}
{\frac{1}{\sigma_{m \mid g}^{2}} - \frac{3\left( y_{m} - \mu_{m \mid g} \right)^{2}}{\sigma_{m \mid g}^{4}},} & {{\text{if}\mspace{6mu}}m = n,\ g = h,} \\
{0,} & \text{otherwise,}
\end{cases}$$

$$\frac{\partial^{2}l_{y_{m} \mid x_{g}}}{\partial\mu_{m \mid g}\,\partial\sigma_{n \mid h}} = \begin{cases}
{-\frac{2\left( y_{m} - \mu_{m \mid g} \right)}{\sigma_{m \mid g}^{3}},} & {{\text{if}\mspace{6mu}}m = n,\ g = h,} \\
{0,} & \text{otherwise.}
\end{cases}$$

Consequently, the Hessian for each conditional parameter has the
following block form:
$${Hess}\left( l_{y_{m} \mid x_{g}} \right)\; = \;\begin{bmatrix}
\frac{\partial^{2}l_{y_{m} \mid x_{g}}}{\partial\mu_{m \mid g}\,\partial\mu_{n \mid h}} & \frac{\partial^{2}l_{y_{m} \mid x_{g}}}{\partial\mu_{m \mid g}\,\partial\sigma_{n \mid h}} \\
\frac{\partial^{2}l_{y_{m} \mid x_{g}}}{\partial\sigma_{n \mid h}\,\partial\mu_{m \mid g}} & \frac{\partial^{2}l_{y_{m} \mid x_{g}}}{\partial\sigma_{m \mid g}\,\partial\sigma_{n \mid h}}
\end{bmatrix}.$$

Notice that each conditional parameter $l_{y_{m} \mid x_{g}}$ has a
Hessian matrix.

## Model for the latent class probabilities

Probabilities of class membership are parameterized with the softmax
transformation:

$$P\left( x_{g} \right) = \frac{\exp\left( \theta_{g} \right)}{\sum\limits_{j}\exp\left( \theta_{j} \right)},$$
where $\theta_{g}$ is the log-scale parameter associated with class $g$.

The jacobian of this transformation is given by

$$J = {diag}(P) - PP^{\top}.$$ Finally, the Hessian for each probability
is

$${Hess}(P\left( x_{g} \right))\; = \; P\left( x_{g} \right)(\left( e_{g} - P \right)\left( e_{g} - P \right)^{\top} - J),$$
where $e_{g}$ is a vector of zeroes with a $1$ in position $g$.

## Model for the conditional probabilities of the multinomial model

Probabilities of conditional responses are parameterized with the
softmax transformation:

$$\pi_{m_{k} \mid g} = \frac{\exp\left( \eta_{m_{k} \mid g} \right)}{\sum\limits_{j}\exp\left( \eta_{m_{j} \mid g} \right)},$$
where $\eta_{m_{k} \mid g}$ is the log-scale parameter associated with
response $k$ to item $m$ in class $g$.

The jacobian of this transformation is given by

$$J_{m \mid g} = {diag}\left( \pi_{m \mid g} \right) - \pi_{m \mid g}\left( \pi_{m \mid g} \right)^{\top}.$$
Finally, the Hessian for each probability is

$${Hess}(\pi_{m_{k} \mid g})\; = \;\pi_{m \mid g}(\left( e_{k} - \pi_{m \mid g} \right)\left( e_{k} - \pi_{m \mid g} \right)^{\top} - J),$$
where $e_{k}$ is a vector of zeroes with a $1$ in position $k$.

## Constant priors

### For latent class probabilities

### For conditional likelihoods

#### Multinomial

For the conditional probabilities modeled with a multinomial likelihood,
we add the following term to the log-likelihood for each class $g$:

$$\lambda_{2} = \sum\limits_{m}\sum\limits_{k}{\widehat{\pi}}_{m_{k}}\sum\limits_{g}\frac{\alpha}{K}\log\left( \pi_{m_{k} \mid g} \right),$$

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
\lambda_{3} & {= \sum\limits_{g}^{K}\left( -0.5\frac{\alpha}{K}\log(\prod\limits_{j}\sigma_{j \mid g}^{2}) - 0.5\frac{\alpha}{K}\sum\limits_{j}\frac{{\widehat{\sigma}}_{j}^{2}}{\sigma_{j \mid g}^{2}} \right)} \\
 & {= \sum\limits_{g}^{K}\left( -\frac{\alpha}{K}\sum\limits_{j}s_{j \mid g} - 0.5\frac{\alpha}{K}\sum\limits_{j}\frac{{\widehat{\sigma}}_{j}^{2}}{\sigma_{j \mid g}^{2}} \right),}
\end{aligned}$$

where ${\widehat{\sigma}}_{j}^{2}$ is the variance of item $j$.

The first-order derivatives are $$\begin{aligned}
{\frac{\partial\lambda_{3}}{\partial\sigma_{m \mid g}}\;} & {= \;\sum\limits_{g}^{K} - \frac{\alpha}{K\sigma_{m \mid g}} + \frac{\alpha}{K\sigma_{m \mid g}^{3}}{\widehat{\sigma}}_{m \mid g}^{2}} \\
 & {= \sum\limits_{g}^{K}\frac{\alpha}{K}(\frac{{\widehat{\sigma}}_{m \mid g}^{2}}{\sigma_{m \mid g}^{3}} - \frac{1}{\sigma_{m \mid g}}).}
\end{aligned}$$

The second-order derivatives are
$$\frac{\partial^{2}\lambda_{3}}{\partial\sigma_{m \mid g}\partial\sigma_{m \mid g}}\; = \;\sum\limits_{g}^{K}\frac{\alpha}{K}(\frac{1}{\sigma_{m \mid g}} - 3\frac{{\widehat{\sigma}}_{m \mid g}^{2}}{\sigma_{m \mid g}^{4}}).$$
