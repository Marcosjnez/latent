# Factor Model (Diagonally weighted least-squares)

## Parametrization

``` math
R = S -\hat{S}(\theta)
```

## Loss

``` math
loss = 0.5 \, \left \lVert W \odot R \right \rVert^2
```

## Gradient

``` math
\Delta_{\hat{S}(\theta)} = -W \odot R
```

## Differential of the gradient

``` math
\partial\Delta_{\hat{S}(\theta)} = W \odot \partial \hat{S}(\theta)
```

## Hessian

``` math
H_{\hat{S}(\theta)} = \operatorname{diag}(\operatorname{vec}(W))
```
