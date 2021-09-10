# Quadrilateral element for linear 2d elasticity 

## Nomenclature 

### Coordinates of the current point 

``` 
[x , x ] 
  1   2
``` 

### Dimensions of the element 

``` 
[h , h ] 
  1   2
``` 

### Numbering of nodes 

``` 
node[1,1] = [x[1] = 0,x[2] = 0] 

node[1,2] = [x[1] = 0,x[2] = h[2]] 

node[2,1] = [x[1] = h[1],x[2] = 0] 

node[2,2] = [x[1] = h[1],x[2] = h[2]] 

``` 

## Shape functions and interpolation of the displacements 

``` 
N[1,1] = (1-x[1]/h[1])*(1-x[2]/h[2]) 

N[1,2] = ((1-x[1]/h[1])*x[2])/h[2] 

N[2,1] = (x[1]*(1-x[2]/h[2]))/h[1] 

N[2,2] = (x[1]*x[2])/(h[1]*h[2]) 

``` 

Testing numbering of vertices... 
OK 

Testing that interpolated displacements are correct at vertices... 
OK 

We get the following expressions for the displacements 

``` 
u[1] = xi[1]*xi[2]*u[[2,2],1]+xi[1]*(1-xi[2])*u[[2,1],1]+(1-xi[1])*xi[2]*u[[1,2],1]+(1-xi[1])*(1-xi[2])*u[[1,1],1] 

u[2] = xi[1]*xi[2]*u[[2,2],2]+xi[1]*(1-xi[2])*u[[2,1],2]+(1-xi[1])*xi[2]*u[[1,2],2]+(1-xi[1])*(1-xi[2])*u[[1,1],2] 

``` 

## Interpolation of the strains 

``` 
B[1,1,1,1,1] = -1/(2*h[1]) 
B[1,1,1,1,2] = 0 
B[1,1,1,2,1] = -1/(2*h[1]) 
B[1,1,1,2,2] = 0 
B[1,1,2,1,1] = 1/(2*h[1]) 
B[1,1,2,1,2] = 0 
B[1,1,2,2,1] = 1/(2*h[1]) 
B[1,1,2,2,2] = 0 
B[1,2,1,1,1] = -1/(4*h[2]) 
B[1,2,1,1,2] = -1/(4*h[1]) 
B[1,2,1,2,1] = 1/(4*h[2]) 
B[1,2,1,2,2] = -1/(4*h[1]) 
B[1,2,2,1,1] = -1/(4*h[2]) 
B[1,2,2,1,2] = 1/(4*h[1]) 
B[1,2,2,2,1] = 1/(4*h[2]) 
B[1,2,2,2,2] = 1/(4*h[1]) 
B[2,1,1,1,1] = -1/(4*h[2]) 
B[2,1,1,1,2] = -1/(4*h[1]) 
B[2,1,1,2,1] = 1/(4*h[2]) 
B[2,1,1,2,2] = -1/(4*h[1]) 
B[2,1,2,1,1] = -1/(4*h[2]) 
B[2,1,2,1,2] = 1/(4*h[1]) 
B[2,1,2,2,1] = 1/(4*h[2]) 
B[2,1,2,2,2] = 1/(4*h[1]) 
B[2,2,1,1,1] = 0 
B[2,2,1,1,2] = -1/(2*h[2]) 
B[2,2,1,2,1] = 0 
B[2,2,1,2,2] = 1/(2*h[2]) 
B[2,2,2,1,1] = 0 
B[2,2,2,1,2] = -1/(2*h[2]) 
B[2,2,2,2,1] = 0 
B[2,2,2,2,2] = 1/(2*h[2]) 
``` 

## Stiffness matrix 

Testing that elastic energy is retrieved from extracted stiffness matrix... 
OK 

We get the following coefficients of the stiffness operator 
``` 

value = -1/4 
K_λ[1,1,1,1,2,2] = value 
K_λ[1,1,1,2,2,2] = value 
K_λ[1,1,2,2,1,1] = value 
K_λ[1,1,2,2,2,1] = value 
K_λ[1,2,1,1,2,2] = value 
K_λ[1,2,1,2,2,2] = value 
K_λ[1,2,2,1,1,1] = value 
K_λ[1,2,2,1,2,1] = value 
K_λ[2,1,1,1,1,2] = value 
K_λ[2,1,1,2,1,2] = value 
K_λ[2,1,2,2,1,1] = value 
K_λ[2,1,2,2,2,1] = value 
K_λ[2,2,1,1,1,2] = value 
K_λ[2,2,1,2,1,2] = value 
K_λ[2,2,2,1,1,1] = value 
K_λ[2,2,2,1,2,1] = value 

value = 1/4 
K_λ[1,1,1,1,1,2] = value 
K_λ[1,1,1,2,1,2] = value 
K_λ[1,1,2,1,1,1] = value 
K_λ[1,1,2,1,2,1] = value 
K_λ[1,2,1,1,1,2] = value 
K_λ[1,2,1,2,1,2] = value 
K_λ[1,2,2,2,1,1] = value 
K_λ[1,2,2,2,2,1] = value 
K_λ[2,1,1,1,2,2] = value 
K_λ[2,1,1,2,2,2] = value 
K_λ[2,1,2,1,1,1] = value 
K_λ[2,1,2,1,2,1] = value 
K_λ[2,2,1,1,2,2] = value 
K_λ[2,2,1,2,2,2] = value 
K_λ[2,2,2,2,1,1] = value 
K_λ[2,2,2,2,2,1] = value 

value = -h[1]/(3*h[2]) 
K_λ[1,1,2,1,2,2] = value 
K_λ[1,2,2,1,1,2] = value 
K_λ[2,1,2,2,2,2] = value 
K_λ[2,2,2,2,1,2] = value 

value = -h[1]/(6*h[2]) 
K_λ[1,1,2,2,2,2] = value 
K_λ[1,2,2,2,1,2] = value 
K_λ[2,1,2,1,2,2] = value 
K_λ[2,2,2,1,1,2] = value 

value = h[1]/(6*h[2]) 
K_λ[1,1,2,2,1,2] = value 
K_λ[1,2,2,2,2,2] = value 
K_λ[2,1,2,1,1,2] = value 
K_λ[2,2,2,1,2,2] = value 

value = h[1]/(3*h[2]) 
K_λ[1,1,2,1,1,2] = value 
K_λ[1,2,2,1,2,2] = value 
K_λ[2,1,2,2,1,2] = value 
K_λ[2,2,2,2,2,2] = value 

value = -h[2]/(3*h[1]) 
K_λ[1,1,1,2,1,1] = value 
K_λ[1,2,1,2,2,1] = value 
K_λ[2,1,1,1,1,1] = value 
K_λ[2,2,1,1,2,1] = value 

value = -h[2]/(6*h[1]) 
K_λ[1,1,1,2,2,1] = value 
K_λ[1,2,1,2,1,1] = value 
K_λ[2,1,1,1,2,1] = value 
K_λ[2,2,1,1,1,1] = value 

value = h[2]/(6*h[1]) 
K_λ[1,1,1,1,2,1] = value 
K_λ[1,2,1,1,1,1] = value 
K_λ[2,1,1,2,2,1] = value 
K_λ[2,2,1,2,1,1] = value 

value = h[2]/(3*h[1]) 
K_λ[1,1,1,1,1,1] = value 
K_λ[1,2,1,1,2,1] = value 
K_λ[2,1,1,2,1,1] = value 
K_λ[2,2,1,2,2,1] = value 

value = -1/4 
K_μ[1,1,1,2,1,2] = value 
K_μ[1,1,1,2,2,2] = value 
K_μ[1,1,2,1,2,1] = value 
K_μ[1,1,2,2,2,1] = value 
K_μ[1,2,1,1,1,2] = value 
K_μ[1,2,1,1,2,2] = value 
K_μ[1,2,2,1,2,1] = value 
K_μ[1,2,2,2,2,1] = value 
K_μ[2,1,1,2,1,2] = value 
K_μ[2,1,1,2,2,2] = value 
K_μ[2,1,2,1,1,1] = value 
K_μ[2,1,2,2,1,1] = value 
K_μ[2,2,1,1,1,2] = value 
K_μ[2,2,1,1,2,2] = value 
K_μ[2,2,2,1,1,1] = value 
K_μ[2,2,2,2,1,1] = value 

value = 1/4 
K_μ[1,1,1,1,1,2] = value 
K_μ[1,1,1,1,2,2] = value 
K_μ[1,1,2,1,1,1] = value 
K_μ[1,1,2,2,1,1] = value 
K_μ[1,2,1,2,1,2] = value 
K_μ[1,2,1,2,2,2] = value 
K_μ[1,2,2,1,1,1] = value 
K_μ[1,2,2,2,1,1] = value 
K_μ[2,1,1,1,1,2] = value 
K_μ[2,1,1,1,2,2] = value 
K_μ[2,1,2,1,2,1] = value 
K_μ[2,1,2,2,2,1] = value 
K_μ[2,2,1,2,1,2] = value 
K_μ[2,2,1,2,2,2] = value 
K_μ[2,2,2,1,2,1] = value 
K_μ[2,2,2,2,2,1] = value 

value = h[1]/(6*h[2])-(2*h[2])/(3*h[1]) 
K_μ[1,1,1,2,1,1] = value 
K_μ[1,2,1,2,2,1] = value 
K_μ[2,1,1,1,1,1] = value 
K_μ[2,2,1,1,2,1] = value 

value = (-h[2]/(3*h[1]))-h[1]/(6*h[2]) 
K_μ[1,1,1,2,2,1] = value 
K_μ[1,2,1,2,1,1] = value 
K_μ[2,1,1,1,2,1] = value 
K_μ[2,2,1,1,1,1] = value 

value = h[1]/(3*h[2])-h[2]/(3*h[1]) 
K_μ[1,1,2,2,1,2] = value 
K_μ[1,2,2,2,2,2] = value 
K_μ[2,1,2,1,1,2] = value 
K_μ[2,2,2,1,2,2] = value 

value = (-h[2]/(6*h[1]))-h[1]/(3*h[2]) 
K_μ[1,1,2,2,2,2] = value 
K_μ[1,2,2,2,1,2] = value 
K_μ[2,1,2,1,2,2] = value 
K_μ[2,2,2,1,1,2] = value 

value = h[2]/(6*h[1])-(2*h[1])/(3*h[2]) 
K_μ[1,1,2,1,2,2] = value 
K_μ[1,2,2,1,1,2] = value 
K_μ[2,1,2,2,2,2] = value 
K_μ[2,2,2,2,1,2] = value 

value = h[2]/(3*h[1])-h[1]/(3*h[2]) 
K_μ[1,1,1,1,2,1] = value 
K_μ[1,2,1,1,1,1] = value 
K_μ[2,1,1,2,2,1] = value 
K_μ[2,2,1,2,1,1] = value 

value = h[2]/(3*h[1])+(2*h[1])/(3*h[2]) 
K_μ[1,1,2,1,1,2] = value 
K_μ[1,2,2,1,2,2] = value 
K_μ[2,1,2,2,1,2] = value 
K_μ[2,2,2,2,2,2] = value 

value = (2*h[2])/(3*h[1])+h[1]/(3*h[2]) 
K_μ[1,1,1,1,1,1] = value 
K_μ[1,2,1,1,2,1] = value 
K_μ[2,1,1,2,1,1] = value 
K_μ[2,2,1,2,2,1] = value 
``` 
