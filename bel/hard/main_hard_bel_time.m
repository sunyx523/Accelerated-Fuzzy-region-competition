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
load('../../star_noisy.mat') %For noisy image input
I = imresize(I,[512 512]);
[m,n] = size(I);

%Initial level set
u0 = zeros(m,n);
for l = 1:n
    for j = 1:m
        if (sqrt((l - n./2).^2 + (j - m./2)^2.*n./m) - n./4) < 0
%        if (sqrt((i - n.*7./8).^2 + (j - m./2)^2.*n./m) - m./16) < 0
            u0(j,l) = 1.0;
        else
            u0(j,l) = 0.0;
        end
    end
end

%parameters
dx = 1/m;          %Grid point distance
tmax = 100000;
N = 2;           %Dimension
g = 1;           %Regularization weight
dumax = 0.001;
demax = 0.1;
method = 1;
%0:gradient descent
%1:1-order
%2:2-order
%3:semi-implicit

tic
for lambda = [500]                %Fidelity weight
for k = [1]                                %times of actual damping coefficient to theortical best damping coefficient
for b = [2]                          %beta in beltrami regularization
    
 %   figure
    E = zeros(tmax, 1);                     %Initial potential energy function
    Ei = 0;
    DUM = zeros(tmax, 1);
    DUMe = zeros(tmax, 1);
    DUS = zeros(tmax, 1);
    DUSe = zeros(tmax, 1);
    DUMi = 0;
    DUSi = 0;
    DUSA = zeros(tmax/10, 1);
    DDUM = zeros(tmax, 1);
    DDUS = zeros(tmax, 1);
    DE = zeros(tmax, 1);
    DEi = 0;
    u = u0;                                %u for level set(primal variable)
    u_ = u;                                %u_ for previous level set
    dum  = 0;
    dum_ = 0;
    ddum = 0;
    ut = 0;
    ut2 = 0;
    dus  = 0;
    dus_ = 0;
    ddus = 0;
    dusa = 0;
    e = 0;
    e_ = 0;
    de = 0;
    ind = 0;
    t1 = clock;
    num = 1;
    t2 = clock;
    a = 2.*pi.*k.*sqrt(b.*g);
    z = 4*N*b*g/dx/dx;
    
%    while  ((de > demax || ind < 10) && ind < tmax/10)
    while (ind < tmax/10)

        %beta = getBeta(u,dx,b);
        beta = b;
        utemp = u;
        etemp = e;
        dumtemp = dum;
        dustemp = dus;
        dusa = dusa + dus;

%%

        DxF = (u(:,[2:end end],:) - u)./dx;
        DyF = (u([2:end end],:,:) - u)./dx;
        grad = DxF.^2 + DyF.^2;
        dom = 1./sqrt(1 + b.^2.*grad);
        domAvg = sum(sum(dom(:)))./m./n;
        a = 2.*pi.*k.*sqrt(b.*g.*domAvg);     %damping coeffience
%        a = 2.*pi.*k.*sqrt(b.*g);
%        a = getDamp(u,dx,b,g,m,n);
        
%%      
        c1 = sum(sum(u.*I))./sum(sum(u));
        c2 = sum(sum((1 - u).*I))./sum(sum(1 - u));
        Bx = beta.*DxF./(sqrt(1 + beta.^2.*(DxF.^2 + DyF.^2)));
        By = beta.*DyF./(sqrt(1 + beta.^2.*(DxF.^2 + DyF.^2)));
        Du = (Bx - Bx(:,[1 1:end - 1],:))./dx + (By - By([1 1:end - 1],:,:))./dx;
        %dE = - g.*getDu(u,dx,beta) + lambda.*((I - c1).^2 - (I - c2).^2);
%         dE = - g.*Du + lambda.*((I - c1).^2 - (I - c2).^2 ...
%             - 2.*u.*(I - c1).*(I.*sum(sum(u)) - sum(sum(u.*I)).*ones(m,n))./sum(sum(u))./sum(sum(u)) ...
%             + 2.*u.*(I - c2).*(-1.*I.*sum(sum(1 - u)) + ones(m,n).*sum(sum((1 - u).*I)))./sum(sum(1 - u))./sum(sum(1 - u))) ;
        dE = - g.*Du + lambda.*((I - c1).^2 - (I - c2).^2 );
        
        [u,dt] = getU(u,u_,I,method,a,z,lambda,dx,beta,dE,g);
        e = g./beta.*sum(sum(sqrt(1 + beta.^2.*(DxF.^2 + DyF.^2)))) + lambda.*sum(sum(u.*(I - c1).^2)) + lambda.*sum(sum((1 - u).*(I - c2).^2));
%        e = 0.5.*sum(sum(ut2(:)))./dt./dt;
     
        %[e,dE] = getDE(u,I,lambda,dx,beta,g);
%%
        u_ = utemp;
        e_ = etemp;
        t2 = clock;
                 
        if(etime(t2,t1) > ind)
            E(int32(ind*10 + 1)) = e;
            DUM(int32(ind*10 + 1)) = dum;
            DUS(int32(ind*10 + 1)) = dus;
            DDUM(int32(ind*10 + 1)) = ddum;
            DDUS(int32(ind*10 + 1)) = ddus;
            DE(int32(ind*10 + 1)) = de;
            if(mod(ind*10, 10) == 0)
                DUSA(int32(ind + 1)) = dusa;
            end
            ind = ind + 0.1;
        end
        
 %      visualize contour
        if(mod(int32(ind*10),10) == 0)
            imshow(uint8(I.*255)); colormap(gray); hold on
            contour(u, [0.5 0.5], 'r');hold off
            title(sprintf('Contour at Level-Set 0 at Time = %d, beta=%d, lambda=%d', ind, b, lambda));
            drawnow;
        end
        
        num = num + 1;
        ut = abs(u(:) - u_(:));
        ut2 = ut.*ut;
        dus = sum(sum(abs(u(:) - u_(:))));
        dum = max(max(abs(u(:) - u_(:))));
        dus_ = dustemp;
        dum_ = dumtemp;
        ddum = abs(dum - dum_);
        ddus = abs(dus - dus_);
        de = abs(e - e_);
        Ei = [Ei e];
        DEi = [DEi de];
        DUMi = [DUMi dum];
        DUSi = [DUSi dus];

    end
    ind = 1;
    %automatically save the energy
    name = strcat('square512_DUM_hard_b =',num2str(b),',t=',num2str(tmax),',lambda=',num2str(lambda),',k=',num2str(k),'m=',num2str(method),'.mat');
    save(name,'DUM');
    name = strcat('square512_DUS_hard_b =',num2str(b),',t=',num2str(tmax),',lambda=',num2str(lambda),',k=',num2str(k),'m=',num2str(method),'.mat');
    save(name,'DUS');    
    name = strcat('square512_DE_hard_b =',num2str(b),',t=',num2str(tmax),',lambda=',num2str(lambda),',k=',num2str(k),'m=',num2str(method),'.mat');
    save(name,'DE');
    name = strcat('square512_E_hard_b =',num2str(b),',t=',num2str(tmax),',lambda=',num2str(lambda),',k=',num2str(k),'m=',num2str(method),'.mat');
    save(name,'E');    
end
end
end
toc
