function result=wedge(w)
% Wedge operator to transform R3 to SO3
result=[ 0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];

end