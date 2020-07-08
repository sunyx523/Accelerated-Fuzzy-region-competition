% main_hard_bel_time.m
% PDE acceleration for TV regularziation of segmentation
% The segmentation scheme we used here is called fuzzy region competition from:
% Mory-Ardon, 2007, "Fuzzy Region Competition: A Convex Two-Phase Segmentation Framework"
%
% getDE(): Calculate the current potential energy and gradient flow
% getU(): Update level set using acceleration scheme
%
% Yuxin Sun
% 1498290038@qq.com
% Georgia Tech
% 2019.9.5
clear all

image = imread('../../cameraman.png');
%image = rgb2gray(image);
image = double(image);
I = (image - min(image(:)))./(max(image(:)) - min(image(:)));
load('../../square_noisy.mat') %For noisy image input
I = imresize(I,[512 512]);
[m,n] = size(I);

%Initial level set
u0 = zeros(m,n);
for i = 1:n
    for j = 1:m
%        if (sqrt((i - n./2).^2 + (j - m./2)^2.*n./m) - n./4) < 0
        if (sqrt((i - n.*7./8).^2 + (j - m./2)^2.*n./m) - m./16) < 0
            u0(j,i) = 1.0;
        else
            u0(j,i) = 0.0;
        end
    end
end

%parameters
dx = 1/m;          %Grid point distance
tmax = 1000000;
N = 2;           %Dimension
g = 1;           %Regularization weight
dumax = 0.001;
demax = 1;
method = 0;
%0:gradient descent
%1:1-order
%2:2-order
%3:semi-implicit
tic;u2 = NonLinObs_PDE_mex(v,ob1,ob2,ui,1e6,tol);toc;


    

