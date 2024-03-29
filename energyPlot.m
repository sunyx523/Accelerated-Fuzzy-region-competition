% energyPlot.m
% display loglog plot for energy functions
% visualize the result for potential energy variation
%
% Yuxin Sun
% 1498290038@qq.com
% Georgia Tech
% 2019.9.5

filename = 'bel/result/7.7square512/DUS/';%file path
energyFolder = fullfile(strcat(filename));
energyInput = dir(fullfile(energyFolder,'*.mat'));
energyNum = length(energyInput);
load(strcat(filename,energyInput(1).name));
En = DUS;
le= '';

for i=1:energyNum;
    name = strcat(filename,energyInput(i).name);
    energy = load(name);
    energy = energy.DUS;
    En = [En energy];
    name = energyInput(i).name;
    name = name(1:length(name) - 4);
    title = char(length(name) + 3);
    temp = '''';
    title(1) = temp;
    title(2:length(name) + 1) = name;
    title(length(name) + 2) = temp;
    title(length(name) + 3) = ',';
    
    le= strcat(le,title);
end
En = En(:,2:energyNum + 1);
le = le(1:length(le) - 1);
figure
loglog(En,'LineWidth',1);
le = str2num(le(:));
legend(le);
xlabel('CPU Time(0.1s)');
ylabel('Potential Energy');
%legend('acceleration(k=5.5)','primaldual','quadratic(h=0.12)','quadratic(h=0.15)','quadratic(h=0.2)','quadratic(h=0.3)')
%legend('Hard-clip','Chambolle-Pock','h=0.5','h=1','h=3','h=5')
%legend('primal dual','gradient descent','PDE acceleration')
%legend('gradient descent','acceleration with optimal damping','acceleration with 5*optimal damping','primal dual')
%legend('soft h=10','soft h=5','soft h=2','soft h=1','hard')
%legend('new optimal','primal dual','old optimal*5','old optimal')
%legend('b=0.5,lambda=1000','b=0.5,lambda=500','b=0.5,lambda=700','b=1,lambda=1000','b=1,lambda=500','b=1,lambda=700','b=2,lambda=1000','b=2,lambda=500','b=2,lambda=700')
%legend('PDE Acceleration','primal dual')