    function [xx,ww]=gaussla(n)
%Calculation of the zeros and weights of Gauss-Laguerre quadrature % Thomas Vallee March 2009
% thomas.vallee@univ-nantes.fr
    format long;
    b(1)=1;
    for i=1:n % Calculation of the binomial coefficients
        b(i+1)=(factorial(n))/(factorial(i)*(factorial(n+1-i))); 
    end
    for i=1:n+1 % The polynomial coefficients 
        poly(i)=((-1)^(n+1-i))*b(i)*(factorial(n)/(factorial(n+1-i)));
    end
    xx=roots(poly);% The polynomial roots
    for i=1:n % Coefficients of the first derivative of the polynomial
        polycd(i)=poly(i)*(n+1-i); 
    end    
    for i=1:n % Evaluation 
        x=xx(i);
        solde=0; 
        for k=1:n
            solde=solde+polycd(k)*(x^(n-k)); 
        end
        ww(i,1)=(factorial(n)^2)/(xx(i)*(solde^2)); 
    end