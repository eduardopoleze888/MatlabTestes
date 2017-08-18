function [delta, t1, t2, grau, cx]  = pol2roots(a,b,c)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Grau 1
    if a == 0 
        grau = 1;
        delta = 0;
        t1 = -c/b;
        t2 = t1;
        cx = 0;
%% Grau 2
    else
        grau = 2;
        delta = b^2 - 4*a*c;
        %% Raizes Complexas
        if delta < 0
            cx = 1;
            t1 = (-b + sqrt(delta))/(2*a);
            t2 = (-b + sqrt(delta))/(2*a);
        %% Raizes Reais Únicas
        elseif delta == 0
            cx = 0;
            t1 = -b/(2*a);
            t2 = t1;
        %% Raizes Reais Distintas
        else
            cx = 0;
            t1 = (-b + sqrt(delta))/(2*a);
            t2 = (-b + sqrt(delta))/(2*a);
        end
    end
end

