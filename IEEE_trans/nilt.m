%************************ NILT?FUNCTION DEFINITION (basic version) **********************%
function [ft,t]=nilt(F,tm);
alfa=0; M=256; P=2;
N=2*M; wyn=2*P+1;
t=linspace(0,tm,M);
NT=2*tm*N/(N-2); omega=2*pi/NT;
c=alfa+25/NT; s=c-i*omega*(0:N+wyn-2);
Fsc=feval(F,s);
ft=fft(Fsc(1:N)); ft=ft(1:M);
for n=N:N+wyn-2
 ft(n-N+2,:)=Fsc(n+1)*exp(-i*n*omega*t);
end
ft1=cumsum(ft); ft2=zeros(wyn-1,M);
for I=1:wyn-2
ft=ft2+1./diff(ft1);
 ft2=ft1(2:wyn-I,:); ft1=ft;
end
ft=ft2+1./diff(ft1); ft=2*real(ft)-Fsc(1);
ft=exp(c*t)/NT.*ft; ft(1)=2*ft(1);
plot(t,ft); xlabel('t'); ylabel('f(t)'); grid on;
%******************************* F1?SUBFUNCTION DEFINITION ******************************%
function f=F1(s) % subfunction F1 is called like this: nilt(�eF1�f,tm);
 f=F(s); % F(s) transform evaluation
%****************************************************************************************************%