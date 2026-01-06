# Factor Model (Maximum likelihood)

## Parametrization

## Loss

$$loss = \frac{n}{2}\,\left( p\log(2\pi) + \log\det(\widehat{S}(\theta)) + \widehat{S}(\theta)^{-1}S \right)$$

## Gradient

$$\Delta_{\widehat{S}{(\theta)}} = \frac{n}{2}\,\widehat{S}(\theta)^{-1}\left( I - S\;\widehat{S}(\theta)^{-1} \right)$$

## Differential of the gradient

$$\partial\widehat{S}(\theta)^{-1} = -\widehat{S}(\theta)^{-1}\,\partial\widehat{S}(\theta)\,\widehat{S}(\theta)^{-1}$$

$$\partial\Delta_{\widehat{S}{(\theta)}} = \frac{n}{2}\,\partial\widehat{S}(\theta)^{-1}\left( I - S\;\widehat{S}(\theta)^{-1} \right) - \widehat{S}(\theta)^{-1}\left( S\;\partial\widehat{S}(\theta)^{-1} \right)$$

## Hessian

$$\begin{aligned}
{H_{\widehat{S}{(\theta)}} = \;\frac{n}{2}\,} & {\widehat{S}(\theta)^{-1}S\;\widehat{S}(\theta)^{-1} \otimes \widehat{S}(\theta)^{1}\; +} \\
 & {\widehat{S}(\theta)^{1} \otimes \widehat{S}(\theta)^{-1}S\;\widehat{S}(\theta)^{-1}\; -} \\
 & {\widehat{S}(\theta)^{-1} \otimes \widehat{S}(\theta)^{-1}}
\end{aligned}$$
