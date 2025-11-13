# The Latent Class model

Article in progress…

## Latent Class Analysis

Sometimes, people belong to different groups (i.e., classes) that are
due to nonobservable characteristics. This fact conditions their
probability of selecting a particular response option when answering an
item. Latent Class Analysis is a statistical model that estimates the
probability that a person belongs to a particular class and the
conditional probabilities of selecting a particular response option
conditioning in the given class.

### The likelihood

Suppose that a sample of people respond to $J$ items and $\mathbf{y}$ is
a vector that contains the scores to each item $j$. Also, let $K$ denote
the number of latent classes and $x_{k}$, the specific class $k$. Then,
the likelihood of this response pattern $\mathbf{y}$, if it was observed
$n$ times in the sample, can be written as

$$\begin{aligned}
l & {= P(\mathbf{y})^{n}} \\
 & {= (\sum\limits_{k = 1}^{K}P\left( x_{k} \right)P\left( \mathbf{y}|x_{k} \right))^{n},}
\end{aligned}$$

Assuming local independence, we can rewrite the conditional
probabilities as

$$P\left( \mathbf{y}|x_{k} \right) = \prod\limits_{j = 1}^{J}P\left( y_{j}|x_{k} \right),$$

where $y_{j}$ denotes the score in item $j$.

With this assumption, the likelihood can be rewritten as
$$l = (\sum\limits_{k = 1}^{K}P\left( x_{k} \right)\prod\limits_{j = 1}^{J}P\left( y_{j}|x_{k} \right))^{n},$$

and the logarithm likelihood becomes
$$ll = n\log(\sum\limits_{k = 1}^{K}P\left( x_{k} \right)\prod\limits_{j = 1}^{J}P\left( y_{j}|x_{k} \right)).$$

### First-order derivatives

The partial derivative of $ll$ with respect to the probability of
belonging to the class $k$ is
$$\frac{\partial ll}{\partial P\left( x_{g} \right)} = \frac{n}{\sum\limits_{k = 1}^{K}P\left( x_{k} \right)\prod\limits_{j = 1}^{J}P\left( y_{j}|x_{k} \right)}\prod\limits_{j = 1}^{J}P\left( y_{j}|x_{g} \right).$$

On the other hand, the partial derivative of $ll$ with respect to the
probability of scoring a particular $y_{j}$ while belonging to the class
$k$ is
$$\frac{\partial ll}{\partial P\left( y_{m}|x_{g} \right)} = \frac{n}{\sum\limits_{k = 1}^{K}P\left( x_{k} \right)\prod\limits_{j = 1}^{J}P\left( y_{j}|x_{k} \right)}P\left( x_{g} \right)\prod\limits_{j \neq m}P\left( y_{j}|x_{g} \right).$$

### Second-order derivatives

The second partial derivative of $ll$ with respect to the probability of
belonging to the class $k$ is

$$\frac{\partial^{2}ll}{\partial P\left( x_{g} \right)\partial P\left( x_{h} \right)} = -\frac{n\prod\limits_{j = 1}^{J}P\left( y_{j}|x_{h} \right)\prod\limits_{j = 1}^{J}P\left( y_{j}|x_{g} \right)}{(\sum\limits_{k = 1}^{K}P\left( x_{k} \right)\prod\limits_{j = 1}^{J}P\left( y_{j}|x_{k} \right))^{2}}$$

$$\frac{\partial ll}{\partial P\left( y_{m}|x_{g} \right)P\left( y_{m}|x_{g} \right)} = -\frac{nP\left( x_{g} \right)\prod\limits_{j \neq m}P\left( y_{j}|x_{g} \right)P\left( x_{g} \right)\prod\limits_{j \neq m}P\left( y_{j}|x_{g} \right)}{(\sum\limits_{k = 1}^{K}P\left( x_{k} \right)\prod\limits_{j = 1}^{J}P\left( y_{j}|x_{k} \right))^{2}}.$$

$$\frac{\partial ll}{\partial P\left( y_{m}|x_{g} \right)P\left( y_{n}|x_{g} \right)} = \frac{n\sum\limits_{k = 1}^{K}P\left( x_{k} \right)\prod\limits_{j = 1}^{J}P\left( y_{j}|x_{k} \right)P\left( x_{g} \right)\prod\limits_{j \neq m,n}P\left( y_{j}|x_{g} \right) - nP\left( x_{g} \right)\prod\limits_{j \neq n}P\left( y_{j}|x_{g} \right)P\left( x_{g} \right)\prod\limits_{j \neq m}P\left( y_{j}|x_{g} \right)}{(\sum\limits_{k = 1}^{K}P\left( x_{k} \right)\prod\limits_{j = 1}^{J}P\left( y_{j}|x_{k} \right))^{2}}.$$

The second partial derivative of $ll$ with respect to the probability of
scoring a particular $y_{m}$ or $y_{n}$ while belonging to the class $g$
or $h$ is
$$\frac{\partial ll}{\partial P\left( y_{m}|x_{g} \right)P\left( y_{n}|x_{h} \right)} = -\frac{nP\left( x_{h} \right)\prod\limits_{j \neq n}P\left( y_{j}|x_{h} \right)P\left( x_{g} \right)\prod\limits_{j \neq m}P\left( y_{j}|x_{g} \right)}{(\sum\limits_{k = 1}^{K}P\left( x_{k} \right)\prod\limits_{j = 1}^{J}P\left( y_{j}|x_{k} \right))^{2}}.$$

The second partial derivative of $ll$ between the probability of
belonging to the class $g$ and the probability of scoring a particular
$y_{m}$ while belonging to the class $h$ is

$$\frac{\partial ll}{\partial P\left( x_{g} \right)\partial P\left( y_{m}|x_{g} \right)} = \frac{n\prod\limits_{j \neq m}P\left( y_{j}|x_{g} \right)\sum\limits_{k = 1}^{K}P\left( x_{k} \right)\prod\limits_{j = 1}^{J}P\left( y_{j}|x_{k} \right) - n\prod\limits_{j = 1}^{J}P\left( y_{j}|x_{g} \right)P\left( x_{g} \right)\prod\limits_{j \neq m}P\left( y_{j}|x_{g} \right)}{(\sum\limits_{k = 1}^{K}P\left( x_{k} \right)\prod\limits_{j = 1}^{J}P\left( y_{j}|x_{k} \right))^{2}}.$$

$$\frac{\partial ll}{\partial P\left( x_{h} \right)\partial P\left( y_{m}|x_{g} \right)} = -\frac{n\prod\limits_{j = 1}^{J}P\left( y_{j}|x_{h} \right)P\left( x_{g} \right)\prod\limits_{j \neq m}P\left( y_{j}|x_{g} \right)}{(\sum\limits_{k = 1}^{K}P\left( x_{k} \right)\prod\limits_{j = 1}^{J}P\left( y_{j}|x_{k} \right))^{2}}.$$

### Model for the conditional probabilities

#### Bernoulli

When $y_{j}$ is a bernoulli random variable, the conditional probability
becomes
$$P\left( y_{j}|x_{k} \right) = \theta_{j}^{y_{j}}\left( 1 - \theta_{j} \right)^{1 - y_{j}},$$
where $\theta_{j}$ is the probability of endorsing item $j$ (i.e.,
$y_{j} = 1$).

Its partial derivative with respect to $\theta_{j}$ is
$$\frac{\partial P\left( y_{j}|x_{k} \right)}{\theta_{j}} = \frac{1}{2y_{j} - 1}.$$

#### Multinomial

$$\frac{\partial P\left( y_{m}|x_{g} \right)}{\partial\theta_{m_{k}|g}} = \frac{n}{l}P\left( x_{g} \right)\prod\limits_{j \neq m}P\left( y_{j}|x_{g} \right).$$

#### Gaussian

### Link function

softmax

### Evaluating the likelihood

trick to prevent undeflow
