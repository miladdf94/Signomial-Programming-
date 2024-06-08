
# Part 1: Relaxations for signomial geometric programming
# Part 2: Sequential ECP/GP

## General
using JuMP
using MathOptInterface
using LinearAlgebra
using XLSX
using SymPy   # to extract monomials from a polynomial

# Function to collect results and save them in an Excel workbook
function collect_results(model, method, instance)
    # Construct a filename based on method and instance
    filename = "results_$(method)_instance_$(instance).xlsx"
    
    # Add a worksheet
    wb = XLSX.Workbook(filename)
    sheet = XLSX.addsheet(wb, "Results")
    
    # Write data to the worksheet
    XLSX.rename(sheet, "Method:", method)
    XLSX.rename(sheet, "Problem Instance:", instance)
    XLSX.rename(sheet, "Termination Status:", termination_status(model))
    XLSX.rename(sheet, "Objective Function Value:", objective_function(model))
    
    # Save the workbook
    XLSX.save(filename) 
end

# List of methods and their corresponding file paths
methods_dict = Dict(
    "ECP"    => "path/to/ECP_relaxation.jl",
    "other methods"   => "path/to/.... other methods ....",
    "GP-Sub" => "path/to/GP-sub.jl")

# List of problem instances
instance_list = 1:7


### monomial estimator function
function mon_est(u1, u2, x_sol, f, x)
    # Calculate a_i
    a1 = u1(x) / f(x)
    a2 = 0 * u2(x) / f(x)

    # Construct f_hat
    f_hat(x_val) = (u1(x_val) / a1)^a1 * (u2(x_val) / a2)^a2

    println("estimated f_hat(x) =", f_hat(x_sol))  # Print the function f_hat
    return f_hat, a1, a2
end



# Part 1- Solve relaxations and compare 
for instance in instance_list
    for (method, method_path) in methods_dict
        # Read data for problem instances
        include(method_path)
        
        # set solver and optimize
        println(method)
        println(model) 
        set_optimizer(model, optimizer_with_attributes(Mosek.Optimizer))
        # collect and save results
        collect_results(model, method, instance)
    end
end

# Part 2: Sequential Algorithm
for iter=1:max_iter  
    for instance in instance_list
        for jj = 1:2 # i=1 indicates sequential ECP, and i=2 indicates sequential GP.
            if jj == 1 # ECP
                # Part 2.1: solve Sequential ECP ALGORITHM with proposed relaxation as initial solution
                include(methods_dict["ECP"]) #  change if needed
                # for iteration 1, it only solves the relaxation - for iteration > 1, add estimators of subproblems
                if iter > 1
                    for i in 1:N
                        @constraint(model, x[i] - s[i] <= exp(sol_xt[i])*(tilde_x[i]-sol_xt[i]+1))
                    end   
                    for k in K_set
                        for j in Cn[k]
                            @constraint(model, gamma[k,j] - gs[k,j] <= exp(sol_gamma_t[k,j])*(tilde_gamma[k,j]-sol_gamma_t[k,j]+1))
                        end
                    end
                end

            elseif jj == 2 # GP
                # initial sol between x_l and x_u
                # Define the number of random seeds
                num_seeds = 10
                # Loop through each variable and generate random numbers using different seeds
                for i in 1:length(x_l)
                    println("Variable ", i)
                    for seed in 1:num_seeds
                        # Set the seed for random number generation
                        Random.seed!(seed)
                        # Generate a random number between x_l and x_u
                        rand_num = x_l[i] + rand() * (x_u[i] - x_l[i])
                        println("Seed ", seed, ": ", rand_num)
                    end
                end
                # For iteration >= 1, add estimators of subproblems
                include(methods_dict["GP-sub"])
            end

        # Solve the optimization problem
        optimize!(model)
        println("#### OPTIMIZATION Process Iteration ", iter, " #####")
        println("Status: ", termination_status(model))
        println("Objective value: ", objective_value(model))
        println("Optimal x values: ", value.(x))
        println("Optimal tilde_x values: ", value.(tilde_x))
        global sol_xt = value.(tilde_x)
        global sol_gamma_t = value.(tilde_gamma)
        println("Is current sol_xt feasible?: ", iffeasible(sol_xt))
    
        # Collect Statistics: Save results for current iteration
        results[iter, 1] = iter
        results[iter, 2] = termination_status(model)
        results[iter, 3] = objective_value(model)
        results[iter, 4:4+N-1] = value.(x)
        results[iter, N+4:2N+3] = value.(tilde_x)
    
        # Check termination conditions
        if iter > 1
            x_current = results[iter, 4:N+3]
            x_past = results[iter-1, 4:N+3]        
            dist = norm(x_current - x_past)
            results[iter,2N+4] = dist        
            if dist <= epsilon
                println("Optimization converged at iteration $iter with distance $dist.")
                break
            end
        end    
end # Sequential END
