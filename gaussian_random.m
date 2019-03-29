clear all;
close all;
clc;

code_length=1:50;

    dots=3+randn([1,50]);


plot (code_length,dots,'*-');
ax = gca;
ax.FontSize=16;
%ax.YLim = [0.01,1];