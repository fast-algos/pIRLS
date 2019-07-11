using LinearAlgebra
using SparseArrays

# The linear solver used at every iteration.
function solve(A, W, b)
    L = A' * W * A
    return L\b
end

# Finds the initial solution, i.e., the l_2-norm minimizer. Note that for
# both cases with and without constraints the solution looks a little
# different.
function InitialSoln(C,d,A,b)
    if findmax(abs.(C))[1]==0
        return solve(A,I,A'*b)
    else
    	inverseCalc = inv(A'*A)
        v = solve(C',inverseCalc,d-C*inverseCalc*A'*b)
        x = solve(A,I,A'*b+C'*v)
        return x
    end
end

# Returns the next step. The case with and without linear constraints are a little different.
function FindDeltaConstraints(C,d,A,b,x,i,p)
    m = size(A)[1]
    R = Diagonal((abs.(A*x-b)).^(p-2) )   
    s = 0.5*i^((p-2)/p)/m^((p-2)/p)
    g = p*R*(A*x-b)
    R1 = R + s*I
    if findmax(abs.(C))[1]==0
        inverseCalc = solve(A, R1, transpose(A)*g)
        quadform = transpose(g)*A*inverseCalc
         if quadform ==0
         	println("s= ",s)
         	println(norm(inverseCalc,2))
         	println("norm gradient  ", norm(g,2))
         end
        Δ = i*inverseCalc/(2*quadform)
    else
        C_aug = vcat(C,g'*A)
        d_aug = vcat(zeros(size(d)[1]),i/2)
        inverseCalc = inv(A'*R1*A)
        v = solve(C_aug',inverseCalc,d_aug)
        Δ = solve(A,R1,C_aug'*v)
    end
    return Δ
end

# A function that calculates the gradient of ||A(x-scale*delta)-b||_p^p.
# Here A,b are as in the input. We use this in the next function to find a
# scale so that given the current solution x and the next step delta, we
# can scale delta so as to make maximum progress.
function GradientScaledObj(scale,p,z,w)
	v = z-scale*w
	y = abs.(v).^(p-2)
	y1 = v.*y
    return (-1)*(w'*y1)
end

# This finds a scaling so that given the current solution x and the next
# step delta, we can scale delta so as to make maximum progress.
function LineSearchObj(A,b,x,p,Δ)
    L = -1
    U = 1
    w = A*Δ
    z = A*x-b
    while GradientScaledObj(U,p,z,w)<0
        L = U
        U = 2*U
    end
    while GradientScaledObj(L,p,z,w)>0
        U = L
        L = 2*L
    end
    @assert (GradientScaledObj(L,p,z,w) < 0)
    @assert (GradientScaledObj(U,p,z,w) > 0)
    while U-L>1e-4
        if (GradientScaledObj((L+U)/2,p,z,w)>0)
            U = (L+U)/2
        else
            L = (L+U)/2
        end
    end
    α = (L+U)/2
    return α
end

# Evaluates the value of the residual problem.
function EvalResidual(A,b,x,p,Δ)
    R = spdiagm(0=>(abs.(A*x-b)).^(p-2))     
    g = p*R*(A*x-b)
    return transpose(g)*A*Δ - 2*p*p*transpose(A*Δ)*R*A*Δ - (p*norm(A*Δ,p))^p
end

# In our algorithm, we had added a padding to the resistances. This
# function checks if we are making enough progress and appropriately
# adjusts the padding.
function Reduce_i(A,b,x,i,p,Δ)
    R = (abs.(A*x-b)).^(p-2)     
    m = size(A)[1]
    s = 0.5*i^((p-2)/p)/m^((p-2)/p)
    SqTerm = transpose(A*Δ)*((R.+s).*(A*Δ))

    # ratio of the square term and the p-norm term at delta.
    k =p*norm(A*Δ,p)*(norm(A*Δ,p)/(2*p*SqTerm))^(1/(p-1)) 
    λ = 16*p

    # the minimum progress that has to be made, or the minimum approximation factor.
    α0 = findmin([1/(16*λ), 1/((16*λ)^(1/(p-1))*k)])[1]  
    γ = EvalResidual(A,b,x,p,α0*Δ)

    # condition that determines whether the padding i has to be reduced.
    if γ - 0.25*α0*i<=0 || SqTerm >= λ*i 				
        return true
    else
        return false
    end
end


# the main algorithm where parameters are:
# ϵ : accuracy we want to achieve
# A,b : the objective we are minimizing is ||Ax-b||_p^p
# p : the norm we want to minimize
# C,d : The linear constraints are Cx = d
# x : Initial solution
# lb : lower bound on the optimum
function pNorm(ϵ,A,b,p,C,d, x, lb)
	# lb is a lower bound on the objective
    # Initial Solution
    current = norm(A*x - b,p)

    println("initial objective = ",current^p)

    # Check if the initial solution is 0. In that case return 0.
    if current^p ==0					                             	
   	println("Norm = 0")
        println("Max gradient entry:", max_gradient_entry)
    	return x
    end
    iteration = 1
    m = size(b)[1]

    # Initial padding. An upper bound on (Initial_Solution - OPT)/16p.
    i = (current^p - lb^p)/(16*p)

    # Termination condition, if this is achieved, we have a (1+ϵ)-approximate solution				                         
    while i > 2*ϵ*current^p/(16*p*(1+ϵ))
        iteration = iteration+1
        println("Iteration count:", iteration) 

        # Find the next step                 
        Δ = FindDeltaConstraints(C,d,A,b,x,i,p)

        #Find the step size that gives the minimum objective given the step Δ.
        α = LineSearchObj(A,b,x,p,Δ)
        need_to_reduce_i = false

		#Check if we have had sufficient progress, if not reduce the padding
        if Reduce_i(A,b,x,i,p,Δ)
            need_to_reduce_i = true
        end

        # Check if our new norm is less than the current norm, and then
        # update x. It is possible to get a value that does not reduce the
        # norm because the line search does not solve to very high
        # precision.
        if norm(A*(x-α*Δ)-b,p) < current
            x = x-α*Δ
            current = norm(A*x-b,p)
            println("Reducing norm : ", current^p)
        else
        	# If we do not reduce the norm, we reduce the padding i.
            need_to_reduce_i = true
        end

        if need_to_reduce_i
        	println("Reducing i")
            i = i/2
        end
        i = min(i, (current^p - lb^p)/(16*p))
    end
    return x,iteration
end


# ϵ : accuracy we want to achieve
# A,b : the objective we are minimizing is ||Ax-b||_p^p
# p : the norm we want to minimize
# C,d : The linear constraints are Cx = d
function pNorm(ϵ,A,b,p,C,d)
    x = InitialSoln(C,d,A,b)
    m = size(b)[1]
    lb = norm(A*x-b, 2)/m^(1/2-1/p)
    return pNorm(ϵ,A,b,p,C,d,x,lb)#,iter_max)#,res)
end
