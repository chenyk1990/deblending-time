%% Reproducible script for Figure 2
%
close all; clc;clear;

%% setting parameters
dt = 4./1000;
tmax = 2.0;
h1 = [-500:20:500];
h2 = [0:20:1000];
tau = [0.5,1.0,1.6];tau=linspace(0.2,1.8,100);
v = [1500,2400,2300];v=linspace(1500,2800,100);
amp = [1., -1.,1];randn('state',2020);amp=randn(size(tau));
f0 = 20;
snr = 200;
L = 9;
seed=2013;

%% make synthetic data
[d1,h,t] = dbt_hevents(dt,f0,tmax,h1,tau,v,amp,snr,L,seed);
[d2,h,t] = dbt_hevents(dt,f0,tmax,h2,tau,v,amp,snr,L,seed);
[nt,nx]=size(d1);

% figure;
% subplot(1,2,1);imagesc(h1,t,d1);
% subplot(1,2,2);imagesc(h2,t,d2);

%% apply dbt_dithering
randn('state',202122);
shift1=floor(0.1*randn(1,size(d1,2))/dt);   % shift of data1 to data2
shift2=-shift1;                             % shift of data2 to data1

d1shift=dbt_dither(d1,shift1);
d2shift=dbt_dither(d2,shift2);

% figure;
% subplot(1,2,1);imagesc(h1,t,d1shift);
% subplot(1,2,2);imagesc(h2,t,d2shift);

%% blend
d1b=d1+d2shift;
d2b=d2+d1shift;
% figure;
% subplot(1,2,1);imagesc(h1,t,d1b);
% subplot(1,2,2);imagesc(h2,t,d2b);

rand('state',2021222324);
[n1,n2]=size(d1);
D1=zeros(nt,nx);
D2=zeros(nt,nx);
% delta0=delta;
delta=shift2;
delta0=delta;
delta0=delta0+round(48*(rand(1,n2)-0.5));
del00=delta0;    %initial guess (fixed)
del=delta;      %ground truth
d1b=d1+dbt_dither(d2,del);
d2b=d2+dbt_dither(d1,-del);
bd=d1b;
labels={'a)','b)','c)','d)','e)'};
figure('units','normalized','Position',[0.2 0.4 0.4, 1],'color','w');
for niter=1:5
    %     niter=20;
    del0=delta0;del0=del0(:);    %initial guess
    del1=del0;      %updated
    
    subplot(5,1,niter);
    plot([1:n2],del,'k','linewidth',2);hold on;
    plot([1:n2],del0,'b','linewidth',2);
    for iter=1:niter %%outer loop (non-linear iterations)
        del0=del1;
        left=dbt_dither(d2,del0);
        left=dbt_deriv(left, 20, 1, 1);
        n1=size(left,1);
        left=dbt_bandpass(left,0.004,0,5,6,6,0,0);
        A=zeros(n1*n2,n2);
        for i2=1:n2
            A(1+(i2-1)*n1:i2*n1,i2)=left(:,i2);
        end
        bd0=d1+dbt_dither(d2,del0);
        right=bd0-bd;
        right=dbt_bandpass(right,0.004,0,5,6,6,0,0);
        dd=inv(A'*A)*(A'*right(:));
        del1=round(del0+1.0*dd);
    end
    plot([1:n2],del1,'r--','linewidth',2);
    xlim([0,70]);ylim([-50,50]);
    ylabel('Time shifts','Fontsize',14,'fontweight','bold');
    xlabel('Shot NO','Fontsize',14,'fontweight','bold');
    tname=sprintf('Non-linear iteration NO #%d, MSE:%g->%g',niter,sum((del(:)-del00(:)).^2)/nx,sum((del(:)-del1(:)).^2)/nx);
    title(tname,'Fontsize',14,'fontweight','bold');
    set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');    
    xticks(linspace(0,60,7));
    legend('True','Initial','Recovered','Location','Best','Fontsize',12);
    text(-8,90,labels{niter},'color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
    %text=strcat('syn3-',num2str(niter),'-non-iterations');
    %title(text);
    %print(gcf,'-dpng','-r300',strcat(text,'.png'));
    fprintf('Niter=%d,SNR=%g->%g,MSE=%g->%g\n',niter,dbt_snr(del(:),del00(:)),dbt_snr(del(:),del1(:)),sum((del(:)-del00(:)).^2)/nx,sum((del(:)-del1(:)).^2)/nx);
    
end
print(gcf,'-depsc','-r300','fig2.eps');




