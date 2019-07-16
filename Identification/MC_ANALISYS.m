% def. temporary vars
Am_min = DATA(1).Am;
Am_max = DATA(1).Am;
Am_mean = DATA(1).Am;
%
AM_min = DATA(1).AM;
AM_max = DATA(1).AM;
AM_mean = DATA(1).AM;
%
lb_min = DATA(1).lb;
lb_max = DATA(1).lb;
lb_mean = DATA(1).lb;
%
ub_min = DATA(1).ub;
ub_max = DATA(1).ub;
ub_mean = DATA(1).ub;
%
tauS_flat_min = DATA(1).tauS_flat;
tauS_flat_max = DATA(1).tauS_flat;
tauS_flat_mean = DATA(1).tauS_flat;
%
tau_x_min = DATA(1).tau_x;
tau_x_max = DATA(1).tau_x;
tau_x_mean = DATA(1).tau_x;
%
e_tau_min = DATA(1).e_tau;
e_tau_max = DATA(1).e_tau;
e_tau_mean = DATA(1).e_tau;
%
norm_e_tau_min = DATA(1).norm_e_tau;
norm_e_tau_max = DATA(1).norm_e_tau;
norm_e_tau_mean = DATA(1).norm_e_tau;
%
x_min = DATA(1).x;
x_max = DATA(1).x;
x_mean = DATA(1).x;
%
e_x_min = DATA(1).e_x;
e_x_max = DATA(1).e_x;
e_x_mean = DATA(1).e_x;
%
x_std = DATA(1).x;
e_x_std = DATA(1).e_x;
Am_std = DATA(1).Am;
AM_std = DATA(1).AM;
%
nMC = temp;
%
for i = 2:nMC
    Am_min = min(Am_min,DATA(i).Am);
    AM_min = min(AM_min,DATA(i).AM);
    lb_min = min(lb_min,DATA(i).lb);
    ub_min = min(ub_min,DATA(i).ub);
    tauS_flat_min = min(tauS_flat_min,DATA(i).tauS_flat);
    tau_x_min = min(tau_x_min,DATA(i).tau_x);
    e_tau_min = min(e_tau_min,DATA(i).e_tau);
    norm_e_tau_min = min(norm_e_tau_min,DATA(i).norm_e_tau);
    x_min = min(x_min,DATA(i).x);
    e_x_min = min(e_x_min,DATA(i).e_x);
    %
    Am_max = max(Am_max,DATA(i).Am);
    AM_max = max(AM_max,DATA(i).AM);
    lb_max = max(lb_max,DATA(i).lb);
    ub_max = max(ub_max,DATA(i).ub);
    tauS_flat_max = max(tauS_flat_max,DATA(i).tauS_flat);
    tau_x_max = max(tau_x_max,DATA(i).tau_x);
    e_tau_max = max(e_tau_max,DATA(i).e_tau);
    norm_e_tau_max = max(norm_e_tau_max,DATA(i).norm_e_tau);
    x_max = max(x_max,DATA(i).x);
    e_x_max = max(e_x_max,DATA(i).e_x);
    %
    Am_mean = Am_mean+DATA(i).Am;
    AM_mean = AM_mean+DATA(i).AM;
    lb_mean = lb_mean+DATA(i).lb;
    ub_mean = ub_mean+DATA(i).ub;
    tauS_flat_mean = tauS_flat_mean + DATA(i).tauS_flat;
    tau_x_mean = tau_x_mean + DATA(i).tau_x;
    e_tau_mean = e_tau_mean + DATA(i).e_tau;
    norm_e_tau_mean = norm_e_tau_mean + DATA(i).norm_e_tau;
    x_mean = x_mean + DATA(i).x;
    e_x_mean = e_x_mean + DATA(i).e_x;
    %
    x_std = [x_std,DATA(i).x];
    e_x_std = [e_x_std,DATA(i).e_x];
    Am_std = [AM_std,DATA(i).Am];
    AM_std = [AM_std,DATA(i).AM];
    %
end

%
Am_mean = Am_mean./nMC;
AM_mean = AM_mean./nMC;
lb_mean = lb_mean./nMC;
ub_mean = ub_mean./nMC;
tauS_flat_mean = tauS_flat_mean./nMC;
tau_x_mean = tau_x_mean./nMC;
e_tau_mean = e_tau_mean./nMC;
norm_e_tau_mean = norm_e_tau_mean./nMC;
x_mean = x_mean./nMC;
e_x_mean = e_x_mean./nMC;

figure(5)
for k = 1:9
    gca = subplot(330+k);
    set(gca, 'XScale', 'log')
    h1 = fill([1:5:length(lb_min),fliplr(1:5:length(lb_min))],[lb_min(k,1:5:end),fliplr(lb_max(k,1:5:end))], 'r');
    set(h1,'facealpha',.5,'EdgeColor','none')
    hold on     
    h2 = fill([1:5:length(lb_min),fliplr(1:5:length(lb_min))],[ub_min(k,1:5:end),fliplr(ub_max(k,1:5:end))], 'b');
    set(h2,'facealpha',.5,'EdgeColor','none')
    semilogx(1:5:length(lb_min),lb_mean(k,1:5:end),'r', 'LineWidth',2)
    semilogx(1:5:length(ub_min),ub_mean(k,1:5:end),'b', 'LineWidth',2)
    semilogx([1,iter],[gammaR_ref(k),gammaR_ref(k)],'--k','LineWidth',1)
    hold off
    box on
    ylabel(['\theta_',num2str(k)])
end

pp = length(tauS_flat_min)-1;
figure(7)
subplot(211)
    h1 = fill([0:pp,fliplr(0:pp)],[tauS_flat_min,fliplr(tauS_flat_max)], 'k');
    set(h1,'facealpha',.5,'EdgeColor','none')
    hold on
    plot(T,tauS_flat_mean,'k--','LineWidth',1)
    %
    h1 = fill([T,fliplr(T)],[tau_x_min,fliplr(tau_x_max)], 'b');
    set(h1,'facealpha',.5,'EdgeColor','none')
    plot(T,tau_x_mean,'b--','LineWidth',1.5)
    hold off
subplot(212)
    h1 = fill([T,fliplr(T)],[norm_e_tau_min',fliplr(norm_e_tau_max')], 'b');
    set(h1,'facealpha',.5,'EdgeColor','none')
    hold on
    plot(T,norm_e_tau_mean','b','LineWidth',2)
    hold off
    
disp('Am/AM:')
[gammaR_ref, Am_mean,std(Am_std,[],2),AM_mean,std(AM_std,[],2)]
    
disp('x_mean:')
[gammaR_ref,x_mean,std(x_std,[],2),e_x_mean,std(e_x_std,[],2)]



