# Hooke material 

The goal of this maxima script is to compute the expression of the stiffness 
matrix (Mandel representation) of a Hooke material in 2 and 3 dimensions. 

**Note:** to save the output of this script as a markdown file, run the 
following command 

``` 
with_stdout("hooke.md", batchload("hooke.mac")); 
``` 


## 2d elasticity (plane strains) 

``` 
    [ 2 mu + lambda_     lambda_       0   ]
    [                                      ]
C = [    lambda_      2 mu + lambda_   0   ] 
    [                                      ]
    [       0               0         2 mu ]
``` 


## 3d elasticity 

``` 
    [ 2 mu + lambda_     lambda_         lambda_       0     0     0   ]
    [                                                                  ]
    [    lambda_      2 mu + lambda_     lambda_       0     0     0   ]
    [                                                                  ]
    [    lambda_         lambda_      2 mu + lambda_   0     0     0   ]
C = [                                                                  ] 
    [       0               0               0         2 mu   0     0   ]
    [                                                                  ]
    [       0               0               0          0    2 mu   0   ]
    [                                                                  ]
    [       0               0               0          0     0    2 mu ]
``` 
