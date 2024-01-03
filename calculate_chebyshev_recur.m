function ineq_chebyshev_recur = calculate_chebyshev_recur(X,U,T,A_mesh,B_mesh,sigma_w,p)
%% This function generates the uncertainty bound using Hoeffding inequality
    T_recur=3:T;
    ineq_chebyshev_recur=ones(size(A_mesh));
    for i=1:size(T_recur,2)
        T_curr=T_recur(i);
        X_curr=X(1:T_curr+1,:);
        U_curr=U(1:T_curr,:);
    
        ineq_chebyshev_curr = calculate_chebyshev(X_curr,U_curr,T_curr,A_mesh,B_mesh,sigma_w,p);
        ineq_chebyshev_recur=ineq_chebyshev_recur.*ineq_chebyshev_curr;
    end

end