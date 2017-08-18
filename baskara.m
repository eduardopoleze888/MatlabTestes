function [delta, t1, t2, grau, cx]  = pol2roots(a,b,c)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Grau 1
    if a == 0 
        grau = 0;
        delta = 0;
        t1 = -c/b;
        cx = 0;
%% Grau 2
    else
        
        delta = b^2 - 4*a*c;
        t1 = (-b + sqrt(delta))/(2*a);
        t2 = (-b + sqrt(delta))/(2*a);
    end
end

