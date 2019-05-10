clear all;
%Data Preparation
filename = 'Gauss3.dat';
delimiterIn = ' ';
headerlinesIn = 120;
A = importdata(filename,delimiterIn,headerlinesIn);
data = A.data;
x0 = data(:, 2);
y0 = data(:, 1);
m0 = 8;
I = eye(m0);