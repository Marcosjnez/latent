# Factor Model

## Transformation

$$\widehat{S}(\theta) = \Lambda\Psi\Lambda^{\top} + \Theta$$

## Update of the gradient

$$\Delta_{\Lambda} = 2\Delta_{\widehat{S}{(\theta)}}\Lambda\Psi$$

$$\Delta_{\Psi} = \Lambda^{\top}\Delta_{\widehat{S}{(\theta)}}\Lambda$$

$$\Delta_{\Theta} = \Delta_{\widehat{S}{(\theta)}}$$

## Update of the differential

$$\partial\Delta_{\widehat{S}{(\theta)}} = \partial\Lambda\Psi\Lambda^{\top} + \Lambda\partial\Psi\Lambda^{\top} + \Lambda\Psi\partial\Lambda^{\top} + \partial\Theta$$

## Update of the differential of the gradient

$$\partial\Delta_{\Lambda} = 2\left( \partial\Delta_{\widehat{S}{(\theta)}}\Lambda\Psi + \Delta_{\widehat{S}{(\theta)}}(\partial\Lambda\Psi + \Lambda\partial\Psi) \right)$$

$$\partial\Delta_{\Psi} = \partial\Lambda^{\top}\Delta_{\widehat{S}{(\theta)}}\Lambda + \Lambda^{\top}\partial\Delta_{\widehat{S}{(\theta)}}\Lambda + \Lambda^{\top}\Delta_{\widehat{S}{(\theta)}}\partial\Lambda$$

$$\partial\Delta_{\Theta} = \partial\Delta_{\widehat{S}{(\theta)}}$$

## Jacobian

$$J_{\Lambda} = (\Lambda\Psi)^{\top} \otimes I_{p}$$

$$J_{\Psi} = \Lambda \otimes \Lambda$$

$$J_{\Psi} = I_{p^{2} \times p^{2}}$$

## Differential of the jacobian towards the gradient

## Constraints

## Constraints’ derivatives
