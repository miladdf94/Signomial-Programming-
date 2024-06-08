
# This provides data for all instances. For further references on where each problem instance is taken from, please refer to the paper.
##### Functions for finding lower and upper bounds
function lower_bounds(Aminus, x_l, x_u)
      result = 1.0
      for i in 1:(N-2)  # consider only the original variables: exclude variables x_n+1 and x_n+2 for finding bounds on gamma
          if A[i] > 0
              result *= x_l[i]^A[i]
          elseif A[i] < 0
              result *= x_u[i]^(-A[i])
          end
      end
      return result
    end
    # This function is for finding products when we have product of x_u[i] for i \in A^+
    # This is used for finding: upper_x_n+1, upper_x_n+2 and upper_gamma  
function upper_bounds(Aplus, x_l, x_u)
      result = 1.0
      for i in 1:(N-2)  # consider only the original variables: exclude variables x_n+1 and x_n+2 for finding bounds on gamma
          if A[i] > 0
              result *= x_u[i]^A[i]
          elseif A[i] < 0
              result *= x_l[i]^(-A[i])
          end
      end
      return result
end
  
#  PROBLEM DATA: instances are not in order.

if p_ecp == 2
      N = 5  # Number of variables (x1, x2, ..., x5): x4 and x5 are auxiliary variables
      M = 4  # Max number of monomials in a constraint
      K = 4  # Number of signomial constraints
      
      #results = Any[max_iter, 2N+5]
      # Sets
      K_set = range(1,K)
      M_set = Dict(1 => 1:4, 2 => 1:2, 3 => 1:2, 4 => 1:3)  # M[k] for k in K - range of monomials for each constraint
      
      # Sets
      Cn = [[], []]  # Initialize C^-_k for k in K
      Cp = [[], []]  # Initialize C^+_k for k in K
      
      #         monomial 1        monomial 2        monomial 3      monomial 4(j_mk+1)   monomial 5 (j_mk+2)
      #         x1 x2 x3 X4    x1 x2 x3 X4     x1 x2 x3 X4      x1 x2 x3 X4      x1 x2 x3 X4
      a = [ [  [1   1  0  0  0],  [1  1  -1  0  0], [0 0  0 -1  0],  [0  0  1  0 0],  [0  0  0 1] ],   # k: signomial constraint 
            [  [1  -1  0  0  0],  [0  0  0  0  0], [0  0  0  0  0],  [0  0  0 0 0],   [0  0  0  0] ], 
            [  [1   1  0  0  0],  [0  0  0  0  0], [0  0  0  0  0],  [0  0  0 0 0],   [0  0  0  0] ], 
            [  [-1  0  0  1  0],  [-1 0  0  0  0], [0  0  0  0  0],  [0  0  0 0 0],   [0  0  0  0] ]       
            ] # j=5 does not exist so let all be zero
      
      
          # j: monomials 
      #     j1       j2        j3     j4     <= 0  
      c = [ 168.0    3651.2   40000   1   ; # k=0: constraint (obj fun)
            1.0425     -1      0       0  ;
            0.00035    -1      0       0  ;
            1.25       41.63   -1      0
            ] # k=1
      
      c_prime = zeros(size(c))
      d = zeros(N)
      d[N-1] = 1
      d[N] = 0
      #d = [1 , 1, 1] 
      # variable bounds
      big = 100
      #       x1     x2    X3     X4    x5    (x4 and x5 are for objective function remodeling) 
      x_l_org = [ 40     40    60     0.1   1  ]    # we should consider 2 variables for x1 as x1p - x1n :: Think about proper indexing
      x_u_org = [ 44     45    70     0.4   10000     ]      
      
      gamma_l = [168*x_l_org[1]*x_l_org[2]   3651.2*x_l_org[1]*x_l_org[2]/x_u_org[3]     40000/x_u_org[4]    0.1     ;  # k=1 (obj fun)
                 1.0425*x_l_org[1]/x_u_org[2]   0     0      0
                 0.00035*x_l_org[1]*x_l_org[2]  0     0      0
                 1.25*x_l_org[4]/x_u_org[1]     0     0      0]
      gamma_u = [168*x_u_org[1]*x_u_org[2]   3651.2*x_u_org[1]*x_u_org[2]/x_l_org[3]     40000/x_l_org[4]    10000     ;  # k=1 (obj fun)
                 1.0425*x_u_org[1]/x_l_org[2]   0     0      0
                 0.00035*x_u_org[1]*x_u_org[2]  0     0      0
                 1.25*x_u_org[4]/x_l_org[1]     0     0      0]
      #
elseif p_ecp == 1
      N = 4  # Number of variables (x1, x2, ..., x5): x4 and x5 are auxiliary variables
      M = 5  # Max number of monomials in a constraint
      K = 2  # Number of signomial constraints
      
      #results = Any[max_iter, 2N+5]
      # Sets
      K_set = range(1,K)
      M_set = Dict(1 => 1:5, 2 => 1:2)  # M[k] for k in K - range of monomials for each constraint
      
      # Sets
      Cn = [[], []]  # Initialize C^-_k for k in K
      Cp = [[], []]  # Initialize C^+_k for k in K
      
      #         monomial 1        monomial 2        monomial 3      monomial 4(j_mk+1)   monomial 5 (j_mk+2)
      #         x1 x2 x3 X4    x1 x2 x3 X4     x1 x2 x3 X4      x1 x2 x3 X4      x1 x2 x3 X4
      a = [ [  [2  0  0  0 ],  [0  2  0  0 ], [1  1  0  0  ],  [0  0  1  0],  [0  0  0 1] ],   # k: signomial constraint 
            [  [1  1  0  0 ],  [0  0  0  0 ], [0  0  0  0  ],  [0  0  0 0],   [0  0  0  0] ] ] # j=5 does not exist so let all be zero
      
      
          # j: monomials 
      #     j1     j2      j3     j4   j5    <= 0  
      c = [ 6.0    4     -2.5     -1   1; # k=0: constraint (obj fun)
            -1     -8      0       0    0] # k=1
      
      c_prime = zeros(size(c))
      d = zeros(N)
      d[N-1] = 1
      d[N] = -1
      #d = [1 , 1, 1] 
      # variable bounds
      big = 100
      #       x1     x2    X3     X4   x5    (x4 and x5 are for objective function remodeling) 
      x_l = [ 1      1     1      20    20  ]    # we should consider 2 variables for x1 as x1p - x1n :: Think about proper indexing
      x_u = [ 100.0  100   100    100   100     ]      
      
      gamma_l = [0   0     0.1   0.1     ;  # k=1 (obj fun)
                 0.1   0     0    0]
      gamma_u = [0    0   100    100;
                 100   0     0    0]


#
elseif p_ecp == 5
      N = 4  # Number of variables (x1, x2, ..., x5): x4 and x5 are auxiliary variables
      M = 6  # Max number of monomials in a constraint
      K = 2  # Number of signomial constraints
      
      #results = Any[max_iter, 2N+5]
      # Sets
      K_set = range(1,K)
      M_set = Dict(1 => 1:6, 2 => 1:4)  # M[k] for k in K - range of monomials for each constraint
      
      # Sets
      Cn = [[], []]  # Initialize C^-_k for k in K
      Cp = [[], []]  # Initialize C^+_k for k in K
      
      #         monomial 1        monomial 2        monomial 3      monomial 4(j_mk+1)   monomial 5 (j_mk+2)
      #         x1 x2 x3 X4    x1 x2 x3 X4     x1 x2 x3 X4      x1 x2 x3 X4      x1 x2 x3 X4
      a = [ [  [1  0  0  0],  [-1  0  0  0], [0  1  0  0 ],  [0 -1  0  0],  [0  0  -1  0], [0 0 0 1] ],   # k: signomial constraint 
            [  [-1 0  0  0],  [0  -1  0  0], [0  0  1  0 ],  [0  0  0 0],   [0  0  0  0], [0 0 0 0] ] ] # j=5 does not exist so let all be zero
      
      
          # j: monomials 
      #     j1     j2      j3     j4   j5    <= 0  
      c = [ 5.0  5000   46.2    72000   144000   -1; # k=0: constraint (obj fun)
            4    32     120       -1    0        0] # k=1
      
      c_prime = zeros(size(c))
      d = zeros(N)
      d[N-1] = 0
      d[N] = 1
      #d = [1 , 1, 1] 
      # variable bounds
      big = 100
      #       x1     x2    X3     X4   x5    (x4 and x5 are for objective function remodeling) 
      x_l = [ 1      1     1      0.1  ]    # we should consider 2 variables for x1 as x1p - x1n :: Think about proper indexing
      x_u = [ 220   220   220    10000     ]      
      
      gamma_l = [0   0     0   0  0   0.1     ;  # k=1 (obj fun)
                 0   0     0   0  0   0]
      gamma_u = [0   0   0  0  0    10000;
                 0   0  0   0  0   0]
      

elseif p_ecp == 7                 
      N = 5  # Number of variables (x1, x2, ..., x5): x4 and x5 are auxiliary variables
      M = 5  # Max number of monomials in a constraint
      K = 2  # Number of signomial constraints
      
      #results = Any[max_iter, 2N+5]
      # Sets
      K_set = range(1,K)
      M_set = Dict(1 => 1:5, 2 => 1:4)  # M[k] for k in K - range of monomials for each constraint
      
      # Sets
      Cn = [[], []]  # Initialize C^-_k for k in K
      Cp = [[], []]  # Initialize C^+_k for k in K
      
      #         monomial 1        monomial 2        monomial 3      monomial 4(j_mk+1)   monomial 5 (j_mk+2)
      #         x1 x2 x3 X4    x1 x2 x3 X4     x1 x2 x3 X4      x1 x2 x3 X4      x1 x2 x3 X4
      a = [ [  [1  -1  0  0 0],  [1  0  0  0 0], [0  -1  0  0 0],  [0  0  0  1 0],  [0  0  0  0 1] ],   # k: signomial constraint 
            [  [1  0  -1  0 0],  [0  0  1  0 0], [1  0  1  0  0],  [0  0  0 0 0],   [0  0  0  0 0]] ] # j=5 does not exist so let all be zero
      
          # j: monomials 
      #     j1     j2      j3     j4   j5    <= 0  
      c = [ 0.5   -1      -5      -1    1   ; # k=0: constraint (obj fun)
            0.01   0.01  0.0005   -1    0   ] # k=1
      
      c_prime = zeros(size(c))
      d = zeros(N)
      d[N-1] = 1
      d[N] = -1
      #d = [1 , 1, 1] 
      # variable bounds
      big = 100
      #       x1     x2    X3     X4   x5    (x4 and x5 are for objective function remodeling) 
      x_l = [ 70      1     0.5    0.1  0.1 ]    # we should consider 2 variables for x1 as x1p - x1n :: Think about proper indexing
      x_u = [ 1500    30     21    1000  1000    ]      
      
      gamma_l = [0   0.1  0.1  0.1  0     ;  # k=1 (obj fun)
                 0   0    0    0    0]
      gamma_u = [0  1000  1000  1000    0;
                 0    0   0  0   0]
                  
elseif p_ecp == 3
      N = 12  # Number of variables (x1, x2, ..., x5): x4 and x5 are auxiliary variables
      M = 7  # Max number of monomials in a constraint
      K = 5  # Number of signomial constraints
      
      #results = Any[max_iter, 2N+5]
      # Sets
      K_set = range(1,K)
      M_set = Dict(1 => 1:7, 2 => 1:3, 3 => 1:4, 4 => 1:4, 5 => 1:4)  # M[k] for k in K - range of monomials for each constraint
      
      # Sets
      Cn = [[], []]  # Initialize C^-_k for k in K
      Cp = [[], []]  # Initialize C^+_k for k in K
      
      #         monomial 1        monomial 2        monomial 3      monomial 4(j_mk+1)   monomial 5 (j_mk+2)
      #         x1 x2 x3 X4    x1 x2 x3 X4     x1 x2 x3 X4      x1 x2 x3 X4      x1 x2 x3 X4
      
      a = [ [  [0.67  0  0  0  0 0  -0.67 0 0 0 0 0],  [0 0.67  0  0  0  0 0  -0.67 0 0 0 0], [0 0  0 0 0 0 0 0 0 0 0 0],  [1  0  0 0 0 0 0 0 0 0 0 0],  [0 1 0  0 0 0 0 0 0 0 0 0], [0  0 0 0 0 0 0 0 0 0 1 0], [0  0 0 0 0 0 0 0 0 0 0 1] ],   # k: signomial constraint 
            [  [0 0 0 0 1 0 1 0 0 0 0 0],  [1 0  0 0 0 0 0 0 0 0 0 0], [0  0 0 0 0 0 0 0 0 0 0 0],  [0  0 0 0 0 0 0 0 0 0 0 0],   [0  0 0 0 0 0 0 0 0 0 0 0], [0  0 0 0 0 0 0 0 0 0 0 0], [0  0 0 0 0 0 0 0 0 0 0 0 ] ], 
            [  [0 0  0  0  0 1 0 1 0 0 0 0],  [1  0  0 0 0 0 0 0 0 0 0 0], [0  1  0  0 0 0 0 0 0 0 0 0],  [0  0 0 0 0 0 0 0 0 0 0 0],   [0  0 0 0 0 0 0 0 0 0 0 0], [0  0 0 0 0 0 0 0 0 0 0 0], [0  0 0 0 0 0 0 0 0 0 0 0] ], 
            [  [0 0  1 0  -1 0 0 0 0 0  0 0],  [0 0  -0.71  0  -1 0 0 0 0 0 0 0], [0  0  -1.3  0  0 0 1 0 0 0 0 0],  [0  0 0 0 0 0 0 0 0 0 0 0],   [0  0 0 0 0 0 0 0 0 0 0 0], [0  0 0 0 0 0 0 0 0 0 0 0],[0  0 0 0 0 0 0 0 0 0 0 0] ],
            [  [0 0 0 1 0 -1 0 0 0 0 0 0], [0 0 0 -0.71 0 -1 0 0 0 0 0 0 ], [0 0 0 -1.3 0 0 0 1 0 0 0 0], [0  0 0 0 0 0 0 0 0 0 0 0],[0  0 0 0 0 0 0 0 0 0 0 0],[0  0 0 0 0 0 0 0 0 0 0 0], [0  0 0 0 0 0 0 0 0 0 0 0] ]       
            ] # j=5 does not exist so let all be zero
      
      
          # j: monomials 
      #     j1      j2      j3     j4   j5  j6  j7  <= 0  
      c = [ 0.4      0.4     10    -1   -1  -1  1; # k=0: constraint (obj fun)
            0.0588   0.1     -1     0   0   0   0 ;
            0.0588   0.1     0.1   -1   0   0   0  ;
            4        2     0.0588  -1   0   0   0  ;
            4        2     0.0588  -1   0   0   0  ;
            ] # k=1
      
      c_prime = zeros(size(c))
      d = zeros(N)
      d[N-1] = 1
      d[N] = -1
      #d = [1 , 1, 1] 
      # variable bounds
      big = 100
      #       x1     x2    X3     X4    x5    (x4 and x5 are for objective function remodeling) 
      x_l_org = 0.1*ones(1,12)
      x_u_org = 10*ones(1,12)
      x_l_org[11] = 0.1      
      x_l_org[11] = 0.1
      x_u_org[11] = 1000      
      x_u_org[11] = 1000
      x_l = x_l_org
      x_u = x_u_org
      
      high = 1000
      high2 = 0.1
      gamma_l = [high2   high2    0  high2   high2  high2  high2    ;  # k=1 (obj fun)
                 high2   high2    0   0      0     0     0    ;
                 high2   high2  high2   0      0     0     0    ;
                 high2   high2  high2   0      0     0     0    ;                 
                 high2   high2  high2   0      0     0     0    ;                 ]
      gamma_u = [high   high    0  high   high  high  high    ;  # k=1 (obj fun)
                 high   high    0   0      0     0     0    ;
                 high   high  high   0      0     0     0    ;
                 high   high  high   0      0     0     0    ;                 
                 high   high  high   0      0     0     0    ;                 ]

end

if p == 1
      #  PROBLEM 1 DATA
      #  PROBLEM 1 DATA
      N_org = 2  # Number of variables (x1, x2, ..., x5)
      K_org = 2  # Number of signomial constraints/terms
      T = 3  # Max number of monomials in a constraint
      rhs = [0 -8]  # right hand side numbers of constraints (the first one is dummy for obj function)
      #results = Any[max_iter, 2N+5]
      # Sets
      T_set = Dict(1 => 1:3, 2 => 1:1)  # M[k] for k in K - range of monomials for each constraint
      #           m1      m2        m3   
      org_a = [ [  [2  0 ], [0  2 ],  [1  1] ],   # k: signomial constraint 
                  [  [1  1 ], [0  0 ],  [0  0] ] ] # j=5 does not exist so let all be zero

      # modify a with additional variables to transfer obje fun into constraints
      #           m1      m2        m3      monomial 4(j_mk+1)   monomial 5 (j_mk+2)
      a = [ [  [2  0 0 0], [0  2 0 0],  [1  1  0 0],  [0 0 -1 0] , [0 0 0 1] ],   # k: signomial constraint 
            [  [1  1 0 0], [0  0 0 0],  [0  0  0 0],  [0 0 0 0] ,  [0 0 0 0] ] ] # j=5 does not exist so let all be zero

      gamma_par = org_a

      # j: monomials 
      #     j1   j2    j3      j4   j5    <= 0  
      c_organized = [ 6    4     2.5 ; # k=0: constraint (obj fun)
                      1   8     0] # k=1
      delta = [ 1    1    -1 ; # k=0: constraint (obj fun)
                -1   1     0] # k=1

      c = [ 6    4     -2.5 -1  0; # k=0: constraint (obj fun)
            -1   8     0   0   0] # k=1

      #       x1    x2   X3    X4   x5    (x4 and x5 are for objective function remodeling) 
      x_l_org = [ 1     1.0  ]    # we should consider 2 variables for x1 as x1p - x1n :: Think about proper indexing
      x_u_org = [ 10.0  10   ]      

      y_l = [0    0]
      y_u = [log(10.0)    log(10)]
      y_mid = (y_l + y_u)/2
              
          
elseif p == 2
    #  PROBLEM 2 DATA
    N_org = 4  # Number of variables (x1, x2, ..., x5)
    K_org = 3  # Max number of monomials in a constraint
    T = 4  # Number of signomial constraints/terms
    rhs = [0 1 1 1]  # right hand side numbers of constraints (the first one is dummy for obj function)
    
    #results = Any[max_iter, 2N+5]
    # Sets
    T_set = Dict(1 => 1:3, 2 => 1:1, 3 => 1:1, 4 => 1:2)  # M[k] for k in K - range of monomials for each constraint
    
    #           m1              m2          m3      monomial 4(j_mk+1)   monomial 5 (j_mk+2)
    org_a = [ [  [1  1  0  0], [1  1 -1 0],  [0  0  0  -1] ],   # k: signomial constraint 
        [  [1 -1  0  0], [0  0  0 0],  [0  0 0 0] ],
        [  [1  1  0  0], [0  0  0 0],  [0  0 0 0] ],  # j=5 does not exist so let all be zero
        [  [-1 0  0 -1], [-1  0  0 0],  [0  0 0 0] ] ]
    gamma_par = org_a
    
        # j: monomials 
    #     j1      j2      j3      j4   j5    <= 0  
    c_organized = [ 168     3651.2  40000    ; # k=0: constraint (obj fun)
        1.0425   0      0        ;
        0.00035  0      0  ;
        1.25     41.63  0] # k=1
        
    
    delta = [ 1   1   1    ; # k=0: constraint (obj fun)
            1   0   0        ;
            1   0   0  ;
            1   1   0   ] # k=1
    
    #       x1    x2   X3    X4   x5    (x4 and x5 are for objective function remodeling) 
    x_l_org = [ 40    40   60  0.1]    # we should consider 2 variables for x1 as x1p - x1n :: Think about proper indexing
    x_u_org = [ 44    45   70   1.4  ]      
    
    y_l = [ log(40)    log(40)   log(60)  log(0.1)]    # 
    y_u = [ log(44)    log(45)   log(70)   log(1.4)  ]      
    y_mid = (y_l + y_u)/2
    
    ### functions
    
elseif p == 3
      N_org = 8  # Number of variables (x1, x2, ..., x5)
      K_org = 5  # Max number of monomials in a constraint
      T = 5  # Number of signomial constraints/terms
      rhs = [0 1 1 1 1]  # right hand side numbers of constraints (the first one is dummy for obj function)
      
      #results = Any[max_iter, 2N+5]
      # Sets
      T_set = Dict(1 => 1:4, 2 => 1:2, 3 => 1:3, 4 => 1:3, 5 => 1:3)  # M[k] for k in K - range of monomials for each constraint
      org_a = [ [  [0.67  0  0  0  0 0  -0.67 0],  [0 0.67  0  0  0  0 0  -0.67], [1  0  0 0 0 0 0 0],  [0 1 0  0 0 0 0 0], [0  0 0 0 0 0 0 0]],   # k: signomial constraint 
            [  [0 0 0 0 1 0 1 0],  [1 0  0 0 0 0 0 0], [0  0 0 0 0 0 0 0],  [0  0 0 0 0 0 0 0 ],   [0  0 0 0 0 0 0 0], ], 
            [  [0 0  0  0  0 1 0 1],  [1  0  0 0 0 0 0 0], [0  1  0  0 0 0 0 0],  [0  0 0 0 0 0 0 0 ],   [0  0 0 0 0 0 0 0 ] ], 
            [  [0 0  1 0  -1 0 0 0],  [0 0  -0.71  0  -1 0 0 0], [0  0  -1.3  0  0 0 1 0],  [0  0 0 0 0 0 0 0],   [0  0 0 0 0 0 0 0] ],
            [  [0 0 0 1 0 -1 0 0 ], [0 0 0 -0.71 0 -1 0 0 ], [0 0 0 -1.3 0 0 0 1], [0  0 0 0 0 0 0 0],[0  0 0 0 0 0 0 0] ]       
            ] # j=5 does not exist so let all be zero
      gamma_par = org_a
      
          # j: monomials 
      c_organized = [ 0.4      0.4   1   1  10; # k=0: constraint (obj fun)
                   0.0588   0.1     1     0   0   ;
                   0.0588   0.1     0.1   1   0    ;
                   4        2     0.0588  1   0    ;
                  4        2     0.0588   1   0 ;
          ] # k=1
        
      delta = [ 1      1   -1   -1  1; # k=0: constraint (obj fun)
          1   1     -1    0   0   ;
          1   1     1    -1   0    ;
          1   1     1    -1   0    ;
          1   1     1    -1   0 ;
             ] # k=1
      
      #       x1    x2   X3    X4   x5    (x4 and x5 are for objective function remodeling) 
      x_l_org = 0.1*ones(1,8)
      x_u_org = 10*ones(1,8)      
      
      y_l = log(0.1)*ones(1,8)    # 
      y_u = log(10)*ones(1,8)  
      y_mid = (y_l + 2*y_u)/2
            ### functions
elseif p == 5
      
      N_org = 3  # Number of variables (x1, x2, ..., x5)
      K_org = 2  # Max number of monomials in a constraint
      T = 5  # Number of signomial constraints/terms
      rhs = [0 1]  # right hand side numbers of constraints (the first one is dummy for obj function)
      
      #results = Any[max_iter, 2N+5]
      # Sets
      T_set = Dict(1 => 1:5, 2 => 1:3)  # M[k] for k in K - range of monomials for each constraint
      
      a = [ [  [1  0  0],  [-1  0  0], [0  1  0],  [0 -1  0],  [0  0  -1] ],   # k: signomial constraint 
            [  [-1 0  0],  [0  -1  0], [0  0  1],  [0  0  0],   [0  0  0]] ] # j=5 does not exist so let all be zero
      
      
          # j: monomials 
      #     j1     j2      j3     j4   j5    <= 0  
      c_organized = [ 5.0  5000   46.2    72000   144000; # k=0: constraint (obj fun)
                       4    32     120       1    0     ] # k=1
      delta = [ 1  1   1    1   1; # k=0: constraint (obj fun)
                1  1   1   -1    0     ] # k=1
      
      gamma_par = org_a
      # j: monomials 
      
      #       x1    x2   X3    X4   x5    (x4 and x5 are for objective function remodeling) 
      x_l_org = ones(1,3)
      x_u_org = 220*ones(1,3)      
      
      y_l = log(1)*ones(1,3)    # 
      y_u = log(220)*ones(1,3)  
      y_mid = (y_l + y_u)/2     
            
end # of loop for problem 3 data assignemnt r



