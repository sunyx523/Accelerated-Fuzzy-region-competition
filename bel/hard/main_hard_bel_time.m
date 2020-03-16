% main_acc_bel.m
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
image = double(image);
I = (image - min(image(:)))./(max(image(:)) - min(image(:)));
load('../../cameraman_noisy.mat') %For noisy image input
[m,n] = size(I);

%Initial level set
u0 = zeros(m,n);
for i = 1:n
    for j = 1:m
        if (sqrt((i - n./2).^2 + (j - m./2)^2.*n./m) - n./4) < 0
            u0(j,i) = 1.0;
        else
            u0(j,i) = 0.0;
        end
    end
end

%parameters
lambda = 0.1;   %Fidelity weight
dx = 1;          %Grid point distance
tmax = 100000;
N = 2;           %Dimension
g = 1;           %Regularization weight
dumax = 0.01;
method = 1;
%0:gradient descent
%1:1-order
%2:2-order
%3:semi-implicit

tic

for k = [1]                          %times of actual damping coefficient to theortical best damping coefficient
for b = [0.5]                          %beta in beltrami regularization
for method = [1]
    
    figure
    E = zeros(tmax,1);                     %Initial potential energy function
    u = u0;                                %u for level set(primal variable)
    u_ = u;                                %u_ for previous level set
    ind = 0;
    t1 = clock;
    du  =100;
    num = 1;
    
    while  (du > dumax && ind < tmax/10)
%    while (ind < tmax/10)

        %beta = getBeta(u,dx,b);
        beta = b;
        utemp = u;
        z = 4*N*beta*g/dx/dx;
%%

         DxF = (u(:,[2:end end],:) - u)./dx;
         DyF = (u([2:end end],:,:) - u)./dx;
         grad = DxF.^2 + DyF.^2;
         dom = sqrt(1 + b.^2.*grad);
         domAvg = sum(sum(dom(:)))./m./n;
         a = 2.*pi.*sqrt(b.*g/m/n).*sqrt(domAvg);
%         a = k.*2.*pi.*sqrt(b.*g);     %damping coeffience
%        a = getDamp(u,dx,b,g,m,n);
        
%%      
        c1 = sum(sum(u.*I))./sum(sum(u));
        c2 = sum(sum((1 - u).*I))./sum(sum(1 - u));
        Bx = beta.*DxF./(sqrt(1 + beta.^2.*(DxF.^2 + DyF.^2)));
        By = beta.*DyF./(sqrt(1 + beta.^2.*(DxF.^2 + DyF.^2)));
        Du = (Bx - Bx(:,[1 1:end-1],:))./dx + (By - By([1 1:end-1],:,:))./dx;
        %dE = - g.*getDu(u,dx,beta) + lambda.*((I - c1).^2 - (I - c2).^2);
        dE = - g.*Du + lambda.*((I - c1).^2 - (I - c2).^2);
        e = g./beta.*sum(sum(sqrt(1 + beta.^2.*(DxF.^2 + DyF.^2)))) + lambda.*sum(sum(u.*(I - c1).^2)) + lambda.*sum(sum((1 - u).*(I - c2).^2));
        
        %[e,dE] = getDE(u,I,lambda,dx,beta,g);
%%
        [u,dt] = getU(u,u_,I,method,a,z,lambda,dx,beta,dE,g);
        u_ = utemp;
        t2 = clock;
        
        if(etime(t2,t1) > ind)
            E(int32(ind*10 + 1)) = e;
            ind = ind + 0.1;
        end
        
 %      visualize contour
%         if(mod(int32(ind*10),10) == 0)
%             imshow(uint8(I.*(max(image(:)) - min(image(:))) + min(image(:)))); colormap(gray); hold on
%             contour(u, [0.5 0.5], 'r');hold off
%             title(sprintf('Contour at Level-Set 0 at Time = %d, beta=%d, lambda=%d', ind, b, lambda));
% %            imshow(u); hold on
%             drawnow;
%         end
        num = num + 1;
        du = max(abs(u(:) - u_(:)));
        
    end
    ind = 1;
    %automatically save the energy
    name = strcat('t_hard_b =',num2str(b),',t=',num2str(tmax),',lambda=',num2str(lambda),',k=',num2str(k),'.mat');
    save(name,'E');
    

end
end
end
toc

    

