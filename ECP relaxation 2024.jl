
# ECP relaxations for signomial geometric programming
# written by: Milad Dehghani Filabadi

using JuMP
using MathOptInterface
using LinearAlgebra

global epsilon = 0.000001
# Read data for problem instances
p = 0
p_ecp = 5
method = "ECP"
path_name = "path.../input_data_instances_all.jl"
include(path_name)

# Populate Cp and Cn to form posynomials based on the sign of c[j, k]
Cn = [[], [], [], [], []]  # Initialize C^-_k for k in K
Cp = [[], [], [], [] , []]  # Initialize C^+_k for k in K
for k in K_set
    copy_M = copy(M_set[k])  # Create a copy of M[k]
    for j in copy_M
        if c[k, j] > 0
            push!(Cp[k], j)
        elseif c[k, j] < 0
            push!(Cn[k], j)
            c_prime[k, j] = -c[k, j]
        end
    end
end

# This function is for finding products when we have product of x_l[i] for i \in A^+
# This is used for finding: lower_x_n+1, lower_x_n+2 and lower_gamma  

# function to add valid inequalities
function add_cuts(model, x, tilde_x, gamma)
  # X cut
  for i in 1:N 
      @constraint(model, x[i] <= (x_u[i] - x_l[i]) / (log(x_u[i]) - log(x_l[i])) * (tilde_x[i] - log(x_l[i])) + x_l[i] ) # 24.c
  end
  for k in K_set
      for j in Cn[k]
          if (gamma_l[k,j] != gamma_u[k,j]) && (gamma_l[k,j] != 0) && (gamma_u[k,j] != 0)
              # GAMMA CUT:
              @constraint(model, gamma[k,j] <= (gamma_u[k,j]-gamma_l[k,j])/(log(gamma_u[k,j])-log(gamma_l[k,j]))*(tilde_gamma[k,j]-log(gamma_l[k,j])) + gamma_l[k,j]) # 24.d
          end
      end
  end
end


##### OPTIMIZATION MODEL
model = Model()
@variable(model, x[1:N])
@variable(model, tilde_x[1:N])
@variable(model, gamma[K_set, 1:M])
@variable(model, tilde_gamma[K_set, 1:M])
@variable(model, lambda[K_set, 1:M])

# Objective function
@objective(model, Min, d'x  )  # 11.a

# Constraints
for k in K_set
@constraint(model, sum(c[k,j] * lambda[k,j] for j in Cp[k]) <= sum(c_prime[k,j] * gamma[k, j] for j in Cn[k])) # 11.b
      for j in Cp[k]
            @constraint(model, [sum(a[k][j][i] * tilde_x[i] for i in 1:N); 1 ; lambda[k, j] ] in MOI.ExponentialCone()  ) # 11.c
      end
      for j in Cn[k]
            @constraint(model,  tilde_gamma[k, j] - sum(a[k][j][i] * tilde_x[i] for i in 1:N) <= 0 )  # 11.d
            @constraint(model, [tilde_gamma[k, j]; 1 ; gamma[k, j]] in MOI.ExponentialCone() ) # 11.g
      end
end
for i in 1:N
      @constraint(model, [tilde_x[i]; 1; x[i]] in MOI.ExponentialCone() ) # 11.e
      @constraint(model, x_l[i] <= x[i] <= x_u[i])  # 24.e
      @constraint(model, log(x_l[i]) - epsilon <= tilde_x[i] <= log(x_u[i]) + epsilon)  # 24.f
end

for k in K_set
      for j in Cn[k]
            # bounds for gamma
            @constraint(model, gamma_l[k,j] <= gamma[k,j] <= gamma_u[k,j]) # 24.g
            if (gamma_l[k,j] != gamma_u[k,j]) && (gamma_l[k,j] != 0) && (gamma_u[k,j] != 0)
                  # bounds for gamma_tilde
                  @constraint(model, log(gamma_l[k,j]) - epsilon <= tilde_gamma[k,j] <= log(gamma_u[k,j]) + epsilon) # 24.h
            end
      end
end

###### CALL CUTS
add_cuts(model, x, tilde_x, gamma)

#### Solve the optimization problem
set_optimizer(model, optimizer_with_attributes(Mosek.Optimizer))
#optimize!(model)           



