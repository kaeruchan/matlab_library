clear all;close all;clc;
%box muller

 counter =0;

 %new function generated
  h_1 = 1.0 * sqrt(-2.0*log(randn(1,1))) * cos(2.0*pi*randn(1,1));
  h_2 = 1.0 * sqrt(-2.0*log(randn(1,1))) * sin(2.0*pi*randn(1,1));
  g_1 = 1.0 * sqrt(-2.0*log(randn(1,1))) * cos(2.0*pi*randn(1,1));
  g_2 = 1.0 * sqrt(-2.0*log(randn(1,1))) * sin(2.0*pi*randn(1,1));
  
  
  
  
  
  %plot(x);
  %axis tight;
