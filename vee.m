function result=vee(R)
    % if R(1,2)==-R(2,1) && R(1,3)==-R(3,1) && R(2,3)==-R(3,2)
        result=[R(3,2);R(1,3);R(2,1)];
    % else
    %     error("Vee: Input is not skew-symmetric");
    % end
end