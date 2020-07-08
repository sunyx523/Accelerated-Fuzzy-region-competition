% main_pd_bel_time.m
% Chambolle-Pock primal-dual for Beltrami regularization
% The segmentation scheme we used here is called fuzzy region competition from:
% Mory-Ardon, 2007, "Fuzzy Region Competition: A Convex Two-Phase Segmentation Framework"
% The convex energy is minimized for instance by using Beltrami version Chambolle-Pock  primal-dual algorithm
% Zosso-Bustin, 2014, "A Primal-Dual Projected Gradient Algorithm for Efficient Beltrami Regularization"
%
% getBeta(): An alternative choice for finding best beta for each
%            iteration, but computational cost will increase.
%
% Yuxin Sun
% 1498290038@qq.com
% Georgia Tech
% 2020.4.29
clear all
clc

image = imread('../cameraman.png');
%image = rgb2gray(image);
image = double(image);
I = (image - min(image(:)))./(max(image(:)) - min(image(:)));
load('../star_noisy.mat') %For noisy image input
I = imresize(I,[512 512]);
[m,n] = size(I);

%Initial level set
u0 = zeros(m,n); 
for i = 1:n
    for j = 1:m
        if (sqrt((i - n./2).^2 + (j - m./2)^2.*n./m) - n./4) < 0
%        if (sqrt((i - n.*7./8).^2 + (j - m./2)^2.*n./m) - m./16) < 0
            u0(j,i) = 1.0;
        else
            u0(j,i) = 0.0;
        end
    end
end

px = zeros(m,n); %Initial projector(dual variable) on x direction
py = zeros(m,n); %Initial projector(dual variable) on y dirextion

%parameters
lambda = 500;      %Fidelity weight
dx = 1/m;             %Grid point distance
tmax = 10000;      %Maximum 2wguiterations
L2 = 8;
g = 1;
dumax = 0.01;
demax = 1;

tic
for lambda = [500]
for b = [2]          %beta in beltrami regularization
    
    figure
    a = 2.*pi.*sqrt(b.*g);
    a = 1;
    sigma = dx*sqrt(a/L2); %sigma parameter from the paper
    tau = dx*sqrt(1/a/L2);
    E = zeros(tmax,1);     %Initla potential energy function
    DE = zeros(tmax,1);
    DUM = zeros(tmax,1);
    DUS = zeros(tmax,1);
    u = u0;                %u for level set(primal variable)
    u_ = u;                %u_ for previous level set
    ind = 0;
    e = 0;
    e_ = 0;
    de = 10;
    t1 = clock;
    dum = 0;
    dum_ = 0;
    ddum = 0;
    dus = 0;
    dus_ = 0;
    ddus = 0;
    num = 1;
 
    while (ind < tmax/10)
%    while  (de > demax && ind < tmax/10)
        beta = b;
        utemp = u;
        etemp = e;
        dumtemp = dum;
        dustemp = dus;
        %beta = getBeta(u,dx,b);   
        c1 = sum(sum(u.*I))/sum(sum(u));
        c2 = sum(sum((1 - u).*I))/sum(sum(1 - u));
        r = lambda.*((I - c1).^2 - (I - c2).^2);

        Dux = (u(:,[2:end end],:) - u)./dx;
        Duy = (u([2:end end],:,:) - u)./dx;
        e = 1./beta.*sum(sum(sqrt(1 + beta.^2.*(Dux.^2 + Duy.^2)))) + lambda.*sum(sum(u.*(I - c1).^2)) + lambda.*sum(sum((1 - u).*(I - c2).^2));
        t2 = clock;
        
        if(etime(t2,t1) > ind)
            E(int32(ind.*10 + 1)) = e;
            DUM(int32(ind*10 + 1)) = dum;
            DUS(int32(ind*10 + 1)) = dus;
            DDUM(int32(ind*10 + 1)) = ddum;
            DDUS(int32(ind*10 + 1)) = ddus;            
            DE(int32(ind*10 + 1)) = de;
            ind = ind + 0.1;
        end
        
%        visualize contour
        if(mod(int32(ind*10),10) == 0)
            imshow(uint8(I.*255)); colormap(gray); hold on
            contour(u, [0.5 0.5], 'r');hold off
            title(sprintf('Primal dual Beltrami at Time = %d, beta=%d, lambda=%d', ind, b, lambda));
%            imshow(u); hold on
            drawnow;
        end
        
        %update dual variable
        Den = real(sqrt(beta.^2 - px.^2 - py.^2));
        Den(Den < 0) = 0;
        Px = px - sigma.*beta.*px + beta.*sigma.*Dux.*Den;
        Py = py - sigma.*beta.*py + beta.*sigma.*Duy.*Den;
        px = beta.*Px./max(beta,sqrt(Px.^2 + Py.^2));
        py = beta.*Py./max(beta,sqrt(Px.^2 + Py.^2));
        
        %update primal variable
        Dpx = (px - px(:,[1 1:end-1],:))./dx;
        Dpy = (py - py([1 1:end-1],:,:))./dx;    
        u = max(min(u - tau.*r + tau.*(Dpx + Dpy),1),0);
        u_ = utemp;
        e_ = etemp;
        
        num = num + 1;
        dum = max(max(abs(u(:) - u_(:))));
        dus = sum(sum(abs(u(:) - u_(:))));
        de = abs(e - e_);
        dus_ = dustemp;
        dum_ = dumtemp;
        ddum = abs(dum - dum_);
        ddus = abs(dus - dus_);        
    end
    
    %automatically save the energy

    name = strcat('square512_DUM_pd_b =',num2str(b),',t=',num2str(tmax),',lambda=',num2str(lambda),'.mat');
    save(name,'DUM');
    name = strcat('square512_DUS_pd_b =',num2str(b),',t=',num2str(tmax),',lambda=',num2str(lambda),'.mat');
    save(name,'DUS');    
    name = strcat('square512_DE_pd_b =',num2str(b),',t=',num2str(tmax),',lambda=',num2str(lambda),'.mat');
    save(name,'DE');
    name = strcat('square512_E_pd_b =',num2str(b),',t=',num2str(tmax),',lambda=',num2str(lambda),'.mat');
    save(name,'E');    
    
end
end
toc

% c = u.*c1 + (1 - u).*c2;
% figure
% imshow(c);
