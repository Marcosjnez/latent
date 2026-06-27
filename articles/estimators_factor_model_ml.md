# Factor Model (Maximum likelihood)

## Parametrization

## Loss

``` math
loss = \frac{n}{2} \, \left( p\log(2\pi) + \log\det\big(\hat{S}(\theta)\big) + \hat{S}(\theta)^{-1} S \right)
```

## Gradient

``` math
\Delta_{\hat{S}(\theta)} = \frac{n}{2} \, \hat{S}(\theta)^{-1} \left( I - S \; \hat{S}(\theta)^{-1} \right)
```

## Differential of the gradient

``` math
\partial \hat{S}(\theta)^{-1} = -\hat{S}(\theta)^{-1} \, \partial \hat{S}(\theta) \, \hat{S}(\theta)^{-1}
```

``` math
\partial\Delta_{\hat{S}(\theta)} = \frac{n}{2} \, \partial \hat{S}(\theta)^{-1} \left( I - S \; \hat{S}(\theta)^{-1} \right) - \hat{S}(\theta)^{-1} \left( S \; \partial \hat{S}(\theta)^{-1} \right)
```

## Hessian

``` math
\begin{aligned}
H_{\hat{S}(\theta)} = \;  \frac{n}{2} \, &\hat{S}(\theta)^{-1} S \; \hat{S}(\theta)^{-1} \otimes \hat{S}(\theta)^{1} \;+ \\
&\hat{S}(\theta)^{1} \otimes \hat{S}(\theta)^{-1} S \; \hat{S}(\theta)^{-1} \; - \\
&\hat{S}(\theta)^{-1} \otimes \hat{S}(\theta)^{-1}
\end{aligned}
```
