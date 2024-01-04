function ineq_chisquare = calculate_chisquare(X,U,T,A_mesh,B_mesh,p)
%% This function generates the uncertainty bound using chisquare inequality, assuming standard normal distribution
    X_ii=X(2:T+1,:);
    X_i=X(1:T,:);
    U_i=U;
    Q=sum(X_ii.^2)+A_mesh.^2*sum(X_i.^2)+B_mesh.^2*sum(U_i.^2)-2*A_mesh*(X_ii'*X_i)-2*B_mesh*(X_ii'*U)+2*A_mesh.*B_mesh*(X_i'*U);
    
    t=-log(1-p);
    ineq_chisquare=Q<T+2*sqrt(T*t)+2*t;

end