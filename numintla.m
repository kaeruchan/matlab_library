function [sol]=numintla(fun,nb) 
[xx,w]=gaussla(nb); % Call the gaussla.m file 
sx=size(xx,1);
sol=0;
for i=1:sx
    x=xx(i); 
    fx=eval(fun); 
    sol=sol+w(i)*fx;
end