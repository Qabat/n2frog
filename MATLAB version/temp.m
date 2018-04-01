clc;
clear;
close all; 

test = dlmread('../../fits/YAG 60 5.7574 0.31967.txt');

plot(test(:,1),test(:,2));
hold on
plot(test(:,1),test(:,3));