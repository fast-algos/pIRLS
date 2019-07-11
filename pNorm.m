%The main function. Want to solve min_x ||Ax-b||_p subject to constraints
%Cx = d to a (1+eps) approximation.
function [final,iteration]= pNorm(eps,A,b,p,C,d)
    x = InitialSoln(C,d,A,b);
    m = size(b);
    lb = norm(A*x-b, 2)/m(1)^(1/2-1/p);
    [final,iteration] = pNorm_with_initial_vector(eps,A,b,p,C,d,x,lb);
    fprintf('Final norm %d \n', norm(final,p)^p);
end

% The main algorithm.
function [final_vec,iteration] = pNorm_with_initial_vector(eps,A,b,p,C,d, x, lb)
    % lb is a lower bound on the objective
    % Initial Solution
    current = norm(A*x - b,p); 
    % Check if the initial solution is 0. In that case return 0.
    if current^p ==0					                             
        println("Norm = 0");
        println("Max gradient entry:", max_gradient_entry);
    	iteration = 1;
        return
    end
    iteration = 1;
    m = size(b);
    
    % Initial padding. An upper bound on (Initial_Solution - OPT)/16p.
    i = (current^p - lb^p)/(16*p)	; 
    
    % Termination condition, if this is achieved, we have a (1+eps)-approximate solution
    while i > 2*eps*current^p/(16*p*(1+eps)) 
        iteration = iteration+1;
        fprintf('Iteration count: %d \n', iteration);
        
        % Find the next step
        delta = FindDeltaConstraints(C,d,A,b,x,i,p);  
        
        % Find the step size that gives the minimum objective given the step delta.
        alpha = LineSearchObj(A,b,x,p,delta);           
        need_to_reduce_i = false;
        
        % Check if we have had sufficient progress, if not reduce the padding
        if Reduce_i(A,b,x,i,p,delta)                    
            need_to_reduce_i = true;
        end
        % Check if our new norm is less than the current norm, and then
        % update x. It is possible to get a value that does not reduce the
        % norm because the line search does not solve to very high
        % precision.
        if norm(A*(x-alpha*delta)-b,p) < current        
            x = x-alpha*delta;
            current = norm(A*x-b,p);
            fprintf('Reducing norm: %d \n', current^p);
        else
            % If we do not reduce the norm, we reduce the padding i.
            need_to_reduce_i = true;                    
        end
        if need_to_reduce_i
            i = i/2;
            fprintf('Reducing i: %d \n', i);
        end
        i = min(i, (current^p - lb^p)/(16*p)); 
    end
    final_vec = x;
end

%The linear solver used at every iteration.
function delta = solve(A, w, b)
    m = size(w);
    W = spdiags(w, zeros(1,1), m(1), m(1));
    L = transpose(A) * W * A;
    delta = L\b;
end

%Finds the initial solution, i.e., the l_2-norm minimizer. Note that for
%both cases with and without constraints the solution looks a little
%different.
function soln = InitialSoln(C,d,A,b)
    m = size(A);
    if max(abs(C))==0
        soln = solve(A,ones(m(1)),transpose(A)*b);
    else
    	inverseCalc = inv(transpose(A)*A);
        v = solve(transpose(C),inverseCalc,d-C*inverseCalc*transpose(A)*b);
        soln = solve(A,eye(m(1)),transpose(A)*b+transpose(C)*v);
    end
end


% Returns the next step.
function delta = FindDeltaConstraints(C,d,A,b,x,i,p)      
    m = size(A);
    s = 0.5*i^((p-2)/p)/m(1)^((p-2)/p);
    r = (abs(A*x-b)).^(p-2);
    g = p* r .*(A*x-b);
    r1 = r + s;
    if max(abs(C))==0
        inverseCalc = solve(A, r1, transpose(A)*g);
        quadform = transpose(g)*A*inverseCalc;
        delta = i*inverseCalc/(2*quadform) ; 
    else
        sizeofd = size(d);
        C_aug = [C;transpose(g)*A];
        d_aug = [zeros(sizeofd(1),1);i/2];
        inverseCalc = inv(transpose(A)*R1*A);
        v = solve(transpose(C_aug),inverseCalc,d_aug);
        delta = solve(A,R1,transpose(C_aug)*v);
    end
end


% A function that calculates the gradient of ||A(x-scale*delta)-b||_p^p.
% Here A,b are as in the input. We use this in the next function to find a
% scale so that given the current solution x and the next step delta, we
% can scale delta so as to make maximum progress.
function obj = GradientScaledObj(scale,p,z,w)
    v = z - scale*w;
    y = abs(v).^(p-2);
    y1 = v .* y;
    obj = -1 * (w' * y1);
end
 
% This finds a scaling so that given the current solution x and the next
% step delta, we can scale delta so as to make maximum progress.
function alpha = LineSearchObj(A,b,x,p,delta)
    L = -3;
    U = 3;
    w = A * delta;
    z = A * x - b;
    while GradientScaledObj(U,p,z,w)<0
        L = U;
        U = 2*U;
    end
    while GradientScaledObj(L,p,z,w)>0
        U = L;
        L = 2*L;
    end
    assert (GradientScaledObj(L,p,z,w) < 0);
    assert (GradientScaledObj(U,p,z,w) > 0);
    while abs(U-L)>1e-1
        if (GradientScaledObj((L+U)/2,p,z,w)>0)
            U = (L+U)/2;
        else
            L = (L+U)/2;
        end
    end
    alpha = (L+U)/2;
end


% Evaluates the value of the residual problem.
function gamma = EvalResidual(A,b,x,p,delta)
    m = size(A);
    R = spdiags((abs(A*x-b)).^(p-2),zeros(1,1),m(1),m(1));        
    g = p*R*(A*x-b);	
    gamma = transpose(g)*A*delta - 2*p*p*transpose(A*delta)*R*A*delta - (p*norm(A*delta,p))^p;
end

% In our algorithm, we had added a padding to the resistances. This
% function checks if we are making enough progress and appropriately
% adjusts the padding.
function [out] = Reduce_i(A,b,x,i,p,delta)    
    m = size(A);
    R = abs(A*x-b).^(p-2);
    s = 0.5*i^((p-2)/p)/m(1)^((p-2)/p)	;
    SqTerm = transpose(A*delta)*((R + s) .* (A*delta));
                
    % ratio of the square term and the p-norm term at delta.
    k =p*norm(A*delta,p)*(norm(A*delta,p)/(2*p*SqTerm))^(1/(p-1)) ;
    lambda = 16*p;
                
    % the minimum progress that has to be made, or the minimum approximation factor.
    alpha0 = min(1/(16*lambda), 1/((16*lambda)^(1/(p-1))*k));
    gamma = EvalResidual(A,b,x,p,alpha0*delta);
                
    % condition that determines whether the padding i has to be reduced.
    if gamma < 0.25*alpha0*i || SqTerm >= lambda*i
        out = true;
    else
        out = false;
    end
end
