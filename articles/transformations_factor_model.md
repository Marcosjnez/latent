# Factor Model

## Transformation

``` math
\hat{S}(\theta) = \Lambda \Psi \Lambda^\top + \Theta
```

## Update of the gradient

``` math
\Delta_{\Lambda} = 2\Delta_{\hat{S}(\theta)} \Lambda \Psi
```

``` math
\Delta_{\Psi} = \Lambda^\top \Delta_{\hat{S}(\theta)} \Lambda
```

``` math
\Delta_{\Theta} = \Delta_{\hat{S}(\theta)}
```

## Update of the differential

``` math
\partial \Delta_{\hat{S}(\theta)} = \partial \Lambda \Psi \Lambda^\top +
\Lambda \partial \Psi \Lambda^\top +
\Lambda \Psi \partial \Lambda^\top +
\partial \Theta
```

## Update of the differential of the gradient

``` math
\partial \Delta_{\Lambda} = 2\left( \partial \Delta_{\hat{S}(\theta)} \Lambda \Psi + \Delta_{\hat{S}(\theta)} \big(\partial \Lambda \Psi + \Lambda \partial\Psi \big) \right)
```

``` math
\partial \Delta_{\Psi} = \partial \Lambda^\top \Delta_{\hat{S}(\theta)} \Lambda +
\Lambda^\top \partial \Delta_{\hat{S}(\theta)} \Lambda +
\Lambda^\top \Delta_{\hat{S}(\theta)} \partial \Lambda
```

``` math
\partial \Delta_{\Theta} = \partial \Delta_{\hat{S}(\theta)}
```

## Jacobian

``` math
J_\Lambda = \big(\Lambda \Psi \big)^\top \otimes I_p
```

``` math
J_\Psi = \Lambda \otimes \Lambda
```

``` math
J_\Psi = I_{p^2\times p^2}
```

## Differential of the jacobian towards the gradient

## Constraints

## Constraints’ derivatives
