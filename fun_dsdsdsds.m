clear all;
close all;
clc;

T=0;
for i=0:1:2
    T = @(x) T(x) + exp(-x);
end
T2 = @(x) exp(-x);

I = sum(chebfun(@(x) T(x)*T2(x),[0, Inf]))