# Factor Model (Diagonally weighted least-squares)

## Parametrization

$$R = S - \widehat{S}(\theta)$$

## Loss

$$loss = 0.5\,\left. \parallel W \odot R\parallel \right.^{2}$$

## Gradient

$$\Delta_{\widehat{S}{(\theta)}} = -W \odot R$$

## Differential of the gradient

$$\partial\Delta_{\widehat{S}{(\theta)}} = W \odot \partial\widehat{S}(\theta)$$

## Hessian

$$H_{\widehat{S}{(\theta)}} = \operatorname{diag}\left( \operatorname{vec}(W) \right)$$
