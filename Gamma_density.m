function [y] = Gamma_density(a,m,sigma)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
y=zeros(size(a));
y(a>=m)=(a(a>=m)-m)./(sigma*(sigma+a(a>=m)-m));
end

