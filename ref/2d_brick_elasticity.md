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
node     = [x  = 0, x  = 0] 
    1, 1     1       2

node     = [x  = 0, x  = h ] 
    1, 2     1       2    2

node     = [x  = h , x  = 0] 
    2, 1     1    1   2

node     = [x  = h , x  = h ] 
    2, 2     1    1   2    2

``` 

## Shape functions and interpolation of the displacements 

``` 
             x        x
              1        2
N     = (1 - --) (1 - --) 
 1, 1        h        h
              1        2

             x
              1
        (1 - --) x
             h    2
              1
N     = ----------- 
 1, 2       h
             2

                x
                 2
        x  (1 - --)
         1      h
                 2
N     = ----------- 
 2, 1       h
             1

        x  x
         1  2
N     = ----- 
 2, 2   h  h
         1  2

``` 

Testing numbering of vertices... 
... OK 

Testing that interpolated displacements are correct at vertices... 
... OK 

We get the following expressions for the displacements 

``` 
u  = xi  xi  u          + xi  (1 - xi ) u          + (1 - xi ) xi  u          + (1 - xi ) (1 - xi ) u          
 1     1   2  [2, 2], 1     1        2   [2, 1], 1          1    2  [1, 2], 1          1         2   [1, 1], 1

u  = xi  xi  u          + xi  (1 - xi ) u          + (1 - xi ) xi  u          + (1 - xi ) (1 - xi ) u          
 2     1   2  [2, 2], 2     1        2   [2, 1], 2          1    2  [1, 2], 2          1         2   [1, 1], 2

``` 

## Interpolation of the strains 

                      1
B                = - ---- 
 1, 1, [1, 1], 1     2 h
                        1
B                = 0 
 1, 1, [1, 1], 2
                      1
B                = - ---- 
 1, 1, [1, 2], 1     2 h
                        1
B                = 0 
 1, 1, [1, 2], 2
                    1
B                = ---- 
 1, 1, [2, 1], 1   2 h
                      1
B                = 0 
 1, 1, [2, 1], 2
                    1
B                = ---- 
 1, 1, [2, 2], 1   2 h
                      1
B                = 0 
 1, 1, [2, 2], 2
                      1
B                = - ---- 
 1, 2, [1, 1], 1     4 h
                        2
                      1
B                = - ---- 
 1, 2, [1, 1], 2     4 h
                        1
                    1
B                = ---- 
 1, 2, [1, 2], 1   4 h
                      2
                      1
B                = - ---- 
 1, 2, [1, 2], 2     4 h
                        1
                      1
B                = - ---- 
 1, 2, [2, 1], 1     4 h
                        2
                    1
B                = ---- 
 1, 2, [2, 1], 2   4 h
                      1
                    1
B                = ---- 
 1, 2, [2, 2], 1   4 h
                      2
                    1
B                = ---- 
 1, 2, [2, 2], 2   4 h
                      1
                      1
B                = - ---- 
 2, 1, [1, 1], 1     4 h
                        2
                      1
B                = - ---- 
 2, 1, [1, 1], 2     4 h
                        1
                    1
B                = ---- 
 2, 1, [1, 2], 1   4 h
                      2
                      1
B                = - ---- 
 2, 1, [1, 2], 2     4 h
                        1
                      1
B                = - ---- 
 2, 1, [2, 1], 1     4 h
                        2
                    1
B                = ---- 
 2, 1, [2, 1], 2   4 h
                      1
                    1
B                = ---- 
 2, 1, [2, 2], 1   4 h
                      2
                    1
B                = ---- 
 2, 1, [2, 2], 2   4 h
                      1
B                = 0 
 2, 2, [1, 1], 1
                      1
B                = - ---- 
 2, 2, [1, 1], 2     2 h
                        2
B                = 0 
 2, 2, [1, 2], 1
                    1
B                = ---- 
 2, 2, [1, 2], 2   2 h
                      2
B                = 0 
 2, 2, [2, 1], 1
                      1
B                = - ---- 
 2, 2, [2, 1], 2     2 h
                        2
B                = 0 
 2, 2, [2, 2], 1
                    1
B                = ---- 
 2, 2, [2, 2], 2   2 h
                      2

## Stiffness matrix 

Testing that elastic energy is retrieved from extracted stiffness matrix... 
... OK 
We get the following coefficients of the stiffness operator 

``` 
                        4 h  mu   2 h  mu   2 h  lambda_
                           2         1         2
K                     = ------- + ------- + ------------ 
 [1, 1], 1, [1, 1], 1    3 h       3 h          3 h
                            1         2            1

                        mu   lambda_
K                     = -- + ------- 
 [1, 1], 1, [1, 1], 2   2       2

                        2 h  mu   2 h  mu   h  lambda_
                           2         1       2
K                     = ------- - ------- + ---------- 
 [1, 1], 1, [1, 2], 1    3 h       3 h         3 h
                            1         2           1

                        mu   lambda_
K                     = -- - ------- 
 [1, 1], 1, [1, 2], 2   2       2

                           4 h  mu    h  mu   2 h  lambda_
                              2        1         2
K                     = (- -------) + ----- - ------------ 
 [1, 1], 1, [2, 1], 1       3 h       3 h         3 h
                               1         2           1

                        lambda_   mu
K                     = ------- - -- 
 [1, 1], 1, [2, 1], 2      2      2

                           2 h  mu    h  mu   h  lambda_
                              2        1       2
K                     = (- -------) - ----- - ---------- 
 [1, 1], 1, [2, 2], 1       3 h       3 h        3 h
                               1         2          1

                           mu    lambda_
K                     = (- --) - ------- 
 [1, 1], 1, [2, 2], 2      2        2

                        mu   lambda_
K                     = -- + ------- 
 [1, 1], 2, [1, 1], 1   2       2

                        2 h  mu   4 h  mu   2 h  lambda_
                           2         1         1
K                     = ------- + ------- + ------------ 
 [1, 1], 2, [1, 1], 2    3 h       3 h          3 h
                            1         2            2

                        lambda_   mu
K                     = ------- - -- 
 [1, 1], 2, [1, 2], 1      2      2

                        h  mu   4 h  mu   2 h  lambda_
                         2         1         1
K                     = ----- - ------- - ------------ 
 [1, 1], 2, [1, 2], 2   3 h      3 h          3 h
                           1        2            2

                        mu   lambda_
K                     = -- - ------- 
 [1, 1], 2, [2, 1], 1   2       2

                           2 h  mu    2 h  mu   h  lambda_
                              2          1       1
K                     = (- -------) + ------- + ---------- 
 [1, 1], 2, [2, 1], 2       3 h        3 h         3 h
                               1          2           2

                           mu    lambda_
K                     = (- --) - ------- 
 [1, 1], 2, [2, 2], 1      2        2

                           h  mu    2 h  mu   h  lambda_
                            2          1       1
K                     = (- -----) - ------- - ---------- 
 [1, 1], 2, [2, 2], 2      3 h       3 h         3 h
                              1         2           2

                        2 h  mu   2 h  mu   h  lambda_
                           2         1       2
K                     = ------- - ------- + ---------- 
 [1, 2], 1, [1, 1], 1    3 h       3 h         3 h
                            1         2           1

                        lambda_   mu
K                     = ------- - -- 
 [1, 2], 1, [1, 1], 2      2      2

                        4 h  mu   2 h  mu   2 h  lambda_
                           2         1         2
K                     = ------- + ------- + ------------ 
 [1, 2], 1, [1, 2], 1    3 h       3 h          3 h
                            1         2            1

                           mu    lambda_
K                     = (- --) - ------- 
 [1, 2], 1, [1, 2], 2      2        2

                           2 h  mu    h  mu   h  lambda_
                              2        1       2
K                     = (- -------) - ----- - ---------- 
 [1, 2], 1, [2, 1], 1       3 h       3 h        3 h
                               1         2          1

                        mu   lambda_
K                     = -- + ------- 
 [1, 2], 1, [2, 1], 2   2       2

                           4 h  mu    h  mu   2 h  lambda_
                              2        1         2
K                     = (- -------) + ----- - ------------ 
 [1, 2], 1, [2, 2], 1       3 h       3 h         3 h
                               1         2           1

                        mu   lambda_
K                     = -- - ------- 
 [1, 2], 1, [2, 2], 2   2       2

                        mu   lambda_
K                     = -- - ------- 
 [1, 2], 2, [1, 1], 1   2       2

                        h  mu   4 h  mu   2 h  lambda_
                         2         1         1
K                     = ----- - ------- - ------------ 
 [1, 2], 2, [1, 1], 2   3 h      3 h          3 h
                           1        2            2

                           mu    lambda_
K                     = (- --) - ------- 
 [1, 2], 2, [1, 2], 1      2        2

                        2 h  mu   4 h  mu   2 h  lambda_
                           2         1         1
K                     = ------- + ------- + ------------ 
 [1, 2], 2, [1, 2], 2    3 h       3 h          3 h
                            1         2            2

                        mu   lambda_
K                     = -- + ------- 
 [1, 2], 2, [2, 1], 1   2       2

                           h  mu    2 h  mu   h  lambda_
                            2          1       1
K                     = (- -----) - ------- - ---------- 
 [1, 2], 2, [2, 1], 2      3 h       3 h         3 h
                              1         2           2

                        lambda_   mu
K                     = ------- - -- 
 [1, 2], 2, [2, 2], 1      2      2

                           2 h  mu    2 h  mu   h  lambda_
                              2          1       1
K                     = (- -------) + ------- + ---------- 
 [1, 2], 2, [2, 2], 2       3 h        3 h         3 h
                               1          2           2

                           4 h  mu    h  mu   2 h  lambda_
                              2        1         2
K                     = (- -------) + ----- - ------------ 
 [2, 1], 1, [1, 1], 1       3 h       3 h         3 h
                               1         2           1

                        mu   lambda_
K                     = -- - ------- 
 [2, 1], 1, [1, 1], 2   2       2

                           2 h  mu    h  mu   h  lambda_
                              2        1       2
K                     = (- -------) - ----- - ---------- 
 [2, 1], 1, [1, 2], 1       3 h       3 h        3 h
                               1         2          1

                        mu   lambda_
K                     = -- + ------- 
 [2, 1], 1, [1, 2], 2   2       2

                        4 h  mu   2 h  mu   2 h  lambda_
                           2         1         2
K                     = ------- + ------- + ------------ 
 [2, 1], 1, [2, 1], 1    3 h       3 h          3 h
                            1         2            1

                           mu    lambda_
K                     = (- --) - ------- 
 [2, 1], 1, [2, 1], 2      2        2

                        2 h  mu   2 h  mu   h  lambda_
                           2         1       2
K                     = ------- - ------- + ---------- 
 [2, 1], 1, [2, 2], 1    3 h       3 h         3 h
                            1         2           1

                        lambda_   mu
K                     = ------- - -- 
 [2, 1], 1, [2, 2], 2      2      2

                        lambda_   mu
K                     = ------- - -- 
 [2, 1], 2, [1, 1], 1      2      2

                           2 h  mu    2 h  mu   h  lambda_
                              2          1       1
K                     = (- -------) + ------- + ---------- 
 [2, 1], 2, [1, 1], 2       3 h        3 h         3 h
                               1          2           2

                        mu   lambda_
K                     = -- + ------- 
 [2, 1], 2, [1, 2], 1   2       2

                           h  mu    2 h  mu   h  lambda_
                            2          1       1
K                     = (- -----) - ------- - ---------- 
 [2, 1], 2, [1, 2], 2      3 h       3 h         3 h
                              1         2           2

                           mu    lambda_
K                     = (- --) - ------- 
 [2, 1], 2, [2, 1], 1      2        2

                        2 h  mu   4 h  mu   2 h  lambda_
                           2         1         1
K                     = ------- + ------- + ------------ 
 [2, 1], 2, [2, 1], 2    3 h       3 h          3 h
                            1         2            2

                        mu   lambda_
K                     = -- - ------- 
 [2, 1], 2, [2, 2], 1   2       2

                        h  mu   4 h  mu   2 h  lambda_
                         2         1         1
K                     = ----- - ------- - ------------ 
 [2, 1], 2, [2, 2], 2   3 h      3 h          3 h
                           1        2            2

                           2 h  mu    h  mu   h  lambda_
                              2        1       2
K                     = (- -------) - ----- - ---------- 
 [2, 2], 1, [1, 1], 1       3 h       3 h        3 h
                               1         2          1

                           mu    lambda_
K                     = (- --) - ------- 
 [2, 2], 1, [1, 1], 2      2        2

                           4 h  mu    h  mu   2 h  lambda_
                              2        1         2
K                     = (- -------) + ----- - ------------ 
 [2, 2], 1, [1, 2], 1       3 h       3 h         3 h
                               1         2           1

                        lambda_   mu
K                     = ------- - -- 
 [2, 2], 1, [1, 2], 2      2      2

                        2 h  mu   2 h  mu   h  lambda_
                           2         1       2
K                     = ------- - ------- + ---------- 
 [2, 2], 1, [2, 1], 1    3 h       3 h         3 h
                            1         2           1

                        mu   lambda_
K                     = -- - ------- 
 [2, 2], 1, [2, 1], 2   2       2

                        4 h  mu   2 h  mu   2 h  lambda_
                           2         1         2
K                     = ------- + ------- + ------------ 
 [2, 2], 1, [2, 2], 1    3 h       3 h          3 h
                            1         2            1

                        mu   lambda_
K                     = -- + ------- 
 [2, 2], 1, [2, 2], 2   2       2

                           mu    lambda_
K                     = (- --) - ------- 
 [2, 2], 2, [1, 1], 1      2        2

                           h  mu    2 h  mu   h  lambda_
                            2          1       1
K                     = (- -----) - ------- - ---------- 
 [2, 2], 2, [1, 1], 2      3 h       3 h         3 h
                              1         2           2

                        mu   lambda_
K                     = -- - ------- 
 [2, 2], 2, [1, 2], 1   2       2

                           2 h  mu    2 h  mu   h  lambda_
                              2          1       1
K                     = (- -------) + ------- + ---------- 
 [2, 2], 2, [1, 2], 2       3 h        3 h         3 h
                               1          2           2

                        lambda_   mu
K                     = ------- - -- 
 [2, 2], 2, [2, 1], 1      2      2

                        h  mu   4 h  mu   2 h  lambda_
                         2         1         1
K                     = ----- - ------- - ------------ 
 [2, 2], 2, [2, 1], 2   3 h      3 h          3 h
                           1        2            2

                        mu   lambda_
K                     = -- + ------- 
 [2, 2], 2, [2, 2], 1   2       2

                        2 h  mu   4 h  mu   2 h  lambda_
                           2         1         1
K                     = ------- + ------- + ------------ 
 [2, 2], 2, [2, 2], 2    3 h       3 h          3 h
                            1         2            2

``` 

## Julia code 

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
K_lambda[1,1,1,1,1,1] = (2*h[2])/(3*h[1]) 
K_lambda[1,1,1,1,1,2] = 1/2 
K_lambda[1,1,1,1,2,1] = h[2]/(3*h[1]) 
K_lambda[1,1,1,1,2,2] = -1/2 
K_lambda[1,1,1,2,1,1] = -(2*h[2])/(3*h[1]) 
K_lambda[1,1,1,2,1,2] = 1/2 
K_lambda[1,1,1,2,2,1] = -h[2]/(3*h[1]) 
K_lambda[1,1,1,2,2,2] = -1/2 
K_lambda[1,1,2,1,1,1] = 1/2 
K_lambda[1,1,2,1,1,2] = (2*h[1])/(3*h[2]) 
K_lambda[1,1,2,1,2,1] = 1/2 
K_lambda[1,1,2,1,2,2] = -(2*h[1])/(3*h[2]) 
K_lambda[1,1,2,2,1,1] = -1/2 
K_lambda[1,1,2,2,1,2] = h[1]/(3*h[2]) 
K_lambda[1,1,2,2,2,1] = -1/2 
K_lambda[1,1,2,2,2,2] = -h[1]/(3*h[2]) 
K_lambda[1,2,1,1,1,1] = h[2]/(3*h[1]) 
K_lambda[1,2,1,1,1,2] = 1/2 
K_lambda[1,2,1,1,2,1] = (2*h[2])/(3*h[1]) 
K_lambda[1,2,1,1,2,2] = -1/2 
K_lambda[1,2,1,2,1,1] = -h[2]/(3*h[1]) 
K_lambda[1,2,1,2,1,2] = 1/2 
K_lambda[1,2,1,2,2,1] = -(2*h[2])/(3*h[1]) 
K_lambda[1,2,1,2,2,2] = -1/2 
K_lambda[1,2,2,1,1,1] = -1/2 
K_lambda[1,2,2,1,1,2] = -(2*h[1])/(3*h[2]) 
K_lambda[1,2,2,1,2,1] = -1/2 
K_lambda[1,2,2,1,2,2] = (2*h[1])/(3*h[2]) 
K_lambda[1,2,2,2,1,1] = 1/2 
K_lambda[1,2,2,2,1,2] = -h[1]/(3*h[2]) 
K_lambda[1,2,2,2,2,1] = 1/2 
K_lambda[1,2,2,2,2,2] = h[1]/(3*h[2]) 
K_lambda[2,1,1,1,1,1] = -(2*h[2])/(3*h[1]) 
K_lambda[2,1,1,1,1,2] = -1/2 
K_lambda[2,1,1,1,2,1] = -h[2]/(3*h[1]) 
K_lambda[2,1,1,1,2,2] = 1/2 
K_lambda[2,1,1,2,1,1] = (2*h[2])/(3*h[1]) 
K_lambda[2,1,1,2,1,2] = -1/2 
K_lambda[2,1,1,2,2,1] = h[2]/(3*h[1]) 
K_lambda[2,1,1,2,2,2] = 1/2 
K_lambda[2,1,2,1,1,1] = 1/2 
K_lambda[2,1,2,1,1,2] = h[1]/(3*h[2]) 
K_lambda[2,1,2,1,2,1] = 1/2 
K_lambda[2,1,2,1,2,2] = -h[1]/(3*h[2]) 
K_lambda[2,1,2,2,1,1] = -1/2 
K_lambda[2,1,2,2,1,2] = (2*h[1])/(3*h[2]) 
K_lambda[2,1,2,2,2,1] = -1/2 
K_lambda[2,1,2,2,2,2] = -(2*h[1])/(3*h[2]) 
K_lambda[2,2,1,1,1,1] = -h[2]/(3*h[1]) 
K_lambda[2,2,1,1,1,2] = -1/2 
K_lambda[2,2,1,1,2,1] = -(2*h[2])/(3*h[1]) 
K_lambda[2,2,1,1,2,2] = 1/2 
K_lambda[2,2,1,2,1,1] = h[2]/(3*h[1]) 
K_lambda[2,2,1,2,1,2] = -1/2 
K_lambda[2,2,1,2,2,1] = (2*h[2])/(3*h[1]) 
K_lambda[2,2,1,2,2,2] = 1/2 
K_lambda[2,2,2,1,1,1] = -1/2 
K_lambda[2,2,2,1,1,2] = -h[1]/(3*h[2]) 
K_lambda[2,2,2,1,2,1] = -1/2 
K_lambda[2,2,2,1,2,2] = h[1]/(3*h[2]) 
K_lambda[2,2,2,2,1,1] = 1/2 
K_lambda[2,2,2,2,1,2] = -(2*h[1])/(3*h[2]) 
K_lambda[2,2,2,2,2,1] = 1/2 
K_lambda[2,2,2,2,2,2] = (2*h[1])/(3*h[2]) 
K_mu[1,1,1,1,1,1] = (4*h[2])/(3*h[1])+(2*h[1])/(3*h[2]) 
K_mu[1,1,1,1,1,2] = 1/2 
K_mu[1,1,1,1,2,1] = (2*h[2])/(3*h[1])-(2*h[1])/(3*h[2]) 
K_mu[1,1,1,1,2,2] = 1/2 
K_mu[1,1,1,2,1,1] = h[1]/(3*h[2])-(4*h[2])/(3*h[1]) 
K_mu[1,1,1,2,1,2] = -1/2 
K_mu[1,1,1,2,2,1] = (-(2*h[2])/(3*h[1]))-h[1]/(3*h[2]) 
K_mu[1,1,1,2,2,2] = -1/2 
K_mu[1,1,2,1,1,1] = 1/2 
K_mu[1,1,2,1,1,2] = (2*h[2])/(3*h[1])+(4*h[1])/(3*h[2]) 
K_mu[1,1,2,1,2,1] = -1/2 
K_mu[1,1,2,1,2,2] = h[2]/(3*h[1])-(4*h[1])/(3*h[2]) 
K_mu[1,1,2,2,1,1] = 1/2 
K_mu[1,1,2,2,1,2] = (2*h[1])/(3*h[2])-(2*h[2])/(3*h[1]) 
K_mu[1,1,2,2,2,1] = -1/2 
K_mu[1,1,2,2,2,2] = (-h[2]/(3*h[1]))-(2*h[1])/(3*h[2]) 
K_mu[1,2,1,1,1,1] = (2*h[2])/(3*h[1])-(2*h[1])/(3*h[2]) 
K_mu[1,2,1,1,1,2] = -1/2 
K_mu[1,2,1,1,2,1] = (4*h[2])/(3*h[1])+(2*h[1])/(3*h[2]) 
K_mu[1,2,1,1,2,2] = -1/2 
K_mu[1,2,1,2,1,1] = (-(2*h[2])/(3*h[1]))-h[1]/(3*h[2]) 
K_mu[1,2,1,2,1,2] = 1/2 
K_mu[1,2,1,2,2,1] = h[1]/(3*h[2])-(4*h[2])/(3*h[1]) 
K_mu[1,2,1,2,2,2] = 1/2 
K_mu[1,2,2,1,1,1] = 1/2 
K_mu[1,2,2,1,1,2] = h[2]/(3*h[1])-(4*h[1])/(3*h[2]) 
K_mu[1,2,2,1,2,1] = -1/2 
K_mu[1,2,2,1,2,2] = (2*h[2])/(3*h[1])+(4*h[1])/(3*h[2]) 
K_mu[1,2,2,2,1,1] = 1/2 
K_mu[1,2,2,2,1,2] = (-h[2]/(3*h[1]))-(2*h[1])/(3*h[2]) 
K_mu[1,2,2,2,2,1] = -1/2 
K_mu[1,2,2,2,2,2] = (2*h[1])/(3*h[2])-(2*h[2])/(3*h[1]) 
K_mu[2,1,1,1,1,1] = h[1]/(3*h[2])-(4*h[2])/(3*h[1]) 
K_mu[2,1,1,1,1,2] = 1/2 
K_mu[2,1,1,1,2,1] = (-(2*h[2])/(3*h[1]))-h[1]/(3*h[2]) 
K_mu[2,1,1,1,2,2] = 1/2 
K_mu[2,1,1,2,1,1] = (4*h[2])/(3*h[1])+(2*h[1])/(3*h[2]) 
K_mu[2,1,1,2,1,2] = -1/2 
K_mu[2,1,1,2,2,1] = (2*h[2])/(3*h[1])-(2*h[1])/(3*h[2]) 
K_mu[2,1,1,2,2,2] = -1/2 
K_mu[2,1,2,1,1,1] = -1/2 
K_mu[2,1,2,1,1,2] = (2*h[1])/(3*h[2])-(2*h[2])/(3*h[1]) 
K_mu[2,1,2,1,2,1] = 1/2 
K_mu[2,1,2,1,2,2] = (-h[2]/(3*h[1]))-(2*h[1])/(3*h[2]) 
K_mu[2,1,2,2,1,1] = -1/2 
K_mu[2,1,2,2,1,2] = (2*h[2])/(3*h[1])+(4*h[1])/(3*h[2]) 
K_mu[2,1,2,2,2,1] = 1/2 
K_mu[2,1,2,2,2,2] = h[2]/(3*h[1])-(4*h[1])/(3*h[2]) 
K_mu[2,2,1,1,1,1] = (-(2*h[2])/(3*h[1]))-h[1]/(3*h[2]) 
K_mu[2,2,1,1,1,2] = -1/2 
K_mu[2,2,1,1,2,1] = h[1]/(3*h[2])-(4*h[2])/(3*h[1]) 
K_mu[2,2,1,1,2,2] = -1/2 
K_mu[2,2,1,2,1,1] = (2*h[2])/(3*h[1])-(2*h[1])/(3*h[2]) 
K_mu[2,2,1,2,1,2] = 1/2 
K_mu[2,2,1,2,2,1] = (4*h[2])/(3*h[1])+(2*h[1])/(3*h[2]) 
K_mu[2,2,1,2,2,2] = 1/2 
K_mu[2,2,2,1,1,1] = -1/2 
K_mu[2,2,2,1,1,2] = (-h[2]/(3*h[1]))-(2*h[1])/(3*h[2]) 
K_mu[2,2,2,1,2,1] = 1/2 
K_mu[2,2,2,1,2,2] = (2*h[1])/(3*h[2])-(2*h[2])/(3*h[1]) 
K_mu[2,2,2,2,1,1] = -1/2 
K_mu[2,2,2,2,1,2] = h[2]/(3*h[1])-(4*h[1])/(3*h[2]) 
K_mu[2,2,2,2,2,1] = 1/2 
K_mu[2,2,2,2,2,2] = (2*h[2])/(3*h[1])+(4*h[1])/(3*h[2]) 
``` 