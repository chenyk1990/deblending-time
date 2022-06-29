% script for fig4 (15 minutes?)
%
% This script requires 1) the seistr package (https://github.com/chenyk1990/seistr); 2) the seislet package (currently it is not open-source yet, if you are interested in it, please contact Yangkang Chen, chenyk2016@gmail.com).
close all; clc;clear;
addpath(genpath('./'));
addpath(genpath('~/chenyk/matlibcyk'));
addpath(genpath('~/chenyk/seislet/matfun/src/'));
addpath(genpath('~/seistr'));
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

%% apply dbt_dithering
randn('state',202122);
shift1=floor(0.1*randn(1,size(d1,2))/dt);   % shift of data1 to data2
shift2=-shift1;                             % shift of data2 to data1

d1shift=dbt_dither(d1,shift1);
d2shift=dbt_dither(d2,shift2);

%% blend
d1b=d1+d2shift;
d2b=d2+d1shift;

dip=str_dip2d(d1,5,30,2,0.01, 1, 0.000001,[20,5,1],1);
dip2=str_dip2d(d2,5,30,2,0.01, 1, 0.000001,[20,5,1],1);

%%
% delta0=delta;
rand('state',2021222324);
[n1,n2]=size(d1);
D1=zeros(nt,nx);
D2=zeros(nt,nx);
% delta0=delta;
delta=shift2;
delta0=delta;
% delta0=delta0+round(32*(rand(1,n2)-0.5));
delta0=floor(dbt_meanf(delta,2,1,2));
% delta0(find(abs(delta-delta0)>16))=delta(find(abs(delta-delta0)>16))+round(32*(rand(size(find(abs(delta-delta0)>16)))-0.5));%exclude the outliers
delta0(find(abs(delta-delta0)>12))=delta(find(abs(delta-delta0)>12))+round(24*(rand(size(find(abs(delta-delta0)>12)))-0.5));%exclude the outliers
delta0(43)=delta(43);%guess why? it causes large estimation error in large-offset area
delta0(48)=delta(48);%guess why? it causes large estimation error in large-offset area
del00=delta0;    %initial guess (fixed)

figure;plot(delta);hold on;plot(delta0);legend('Ground truth','Initial','Location','Best');

del=delta;      %ground truth
d1b=d1+dbt_dither(d2,del);
d2b=d2+dbt_dither(d1,-del);
bd=d1b;
snrs=[];
snrs2=[];
mses=[];
data1s=zeros(nt,nx,4);
data2s=zeros(nt,nx,4);
dels=zeros(5,nx);

for iter_out=1:5
    D1=zeros(nt,nx);
    D2=zeros(nt,nx);
    for iter=1:8
        fprintf('\n Iter %d \n',iter);
        D1T=dbt_dither(D1,-delta0);
        D2T=dbt_dither(D2,delta0);
        D1u = D1 + 0.5*(d1b-(D1+D2T));    % updated model
        D2u = D2 + 0.5*(d2b-(D1T+D2));    % updated model
        D1= st_seislet_denoise_2d(D1u,dip,'ps',3,0.1,2,3);
        D2= st_seislet_denoise_2d(D2u,dip2,'ps',3,0.1,2,3);

        snr2(iter)=10*log10(sum(sum(d1.*d1))/sum(sum((d1-D1).*(d1-D1))));
        snr22(iter)=10*log10(sum(sum(d2.*d2))/sum(sum((d2-D2).*(d2-D2))));
%         figure(9);
%         subplot(1,2,1);imagesc(h1,t,D1);
%         subplot(1,2,2);imagesc(h2,t,D2);
        fprintf('iter2=%d,iter1=%d,SNR=%g,SNR2=%g\n',iter_out,iter,snr2(iter),snr22(iter));
    end
    snrs=[snrs,snr2];
    snrs2=[snrs2,snr22];
    data1s(:,:,iter_out)=D1;
    data2s(:,:,iter_out)=D2;
    %% update time
    del0=delta0;del0=del0(:);    %initial guess
    del1=del0;      %updated
    for iter=1:5 %%outer loop (non-linear iterations)
        del0=del1;
        left=dbt_dither(D2,del0);
%           left=dbt_dither(d2,del0);
        left=dbt_deriv(left, 20, 1, 1);
%         left=left(1:10:end,:);
        n1=size(left,1);
        left=dbt_bandpass(left,0.004,0,6,6,6,0,0);
        A=zeros(n1*n2,n2);
        for i2=1:n2
            A(1+(i2-1)*n1:i2*n1,i2)=left(:,i2);
        end
        bd0=D1+dbt_dither(D2,del0);
%           bd0=d1+dbt_dither(d2,del0);
        right=bd0-bd;
%         right=bd0(1:10:end,:)-bd(1:10:end,:);
        right=dbt_bandpass(right,0.004,0,6,6,6,0,0);
        dd=inv(A'*A)*(A'*right(:));
        del1=round(del0+1*dd);
        mmse=sum((del(:)-del1(:)).^2)/nx;
        mses=[mses,mmse];
    end
    dels(iter_out,:)=del1(:)';
    delta0=del1;
fprintf('iter2=%d,SNR=%g->%g,MSE=%g->%g\n',iter_out,yc_snr(del(:),del00(:)),yc_snr(del(:),del1(:)),sum((del(:)-del00(:)).^2)/nx,sum((del(:)-del1(:)).^2)/nx);
end
% figure;plot(snrs);hold on;plot(snrs2);
% figure;plot(mses);

%% maximum db iter2=5,iter1=10,SNR=11.1169,SNR2=10.6301

figure;plot(1:n2,del,'k-','linewidth',2);hold on;plot(1:n2,del00,'b-','linewidth',2);
plot(1:n2,del1,'r--','linewidth',2); 
legend('True','Initial','Recovered','Location','Best');

%% traditional
    D1=zeros(nt,nx);
    D2=zeros(nt,nx);
    for iter=1:8
        fprintf('\n Iter %d \n',iter);
        D1T=dbt_dither(D1,-del00);
        D2T=dbt_dither(D2,del00);
        D1u = D1 + 0.5*(d1b-(D1+D2T));    % updated model
        D2u = D2 + 0.5*(d2b-(D1T+D2));    % updated model
        D1= seislet_denoise_2d(D1u,dip,'ps',3,0.1,2,3);
        D2= seislet_denoise_2d(D2u,dip2,'ps',3,0.1,2,3);
        
        snr2(iter)=10*log10(sum(sum(d1.*d1))/sum(sum((d1-D1).*(d1-D1))));
        snr22(iter)=10*log10(sum(sum(d2.*d2))/sum(sum((d2-D2).*(d2-D2))));
        figure(9);
        subplot(1,2,1);imagesc(h1,t,D1);
        subplot(1,2,2);imagesc(h2,t,D2);
        fprintf('iter2=%d,iter1=%d,SNR=%g,SNR2=%g\n',iter_out,iter,snr2(iter),snr22(iter));
    end
figure;plot(snr2);

figure('units','normalized','Position',[0.2 0.4 0.6, 0.6],'color','w');
subplot(2,1,1);
plot(1:n2,del,'k-','linewidth',2);hold on;plot(1:n2,del00,'b-','linewidth',2);
plot(1:n2,del1,'r--','linewidth',2); 
xlim([0,52]);
xlabel('Shot NO #','Fontsize',12,'fontweight','bold');
ylabel('Time shifts','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
legend('True','Initial','Recovered','Location','Best');
text(-3,65,'a)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','left');

subplot(2,2,3);
plot([1:8*5],snrs,'r','linewidth',2);hold on;plot([1:8*5],snrs2,'b','linewidth',2);
plot([1:8],snr2,'g--','linewidth',2);
plot([1:8],snr22,'c--','linewidth',2);
xlabel('Iteration NO #','Fontsize',12,'fontweight','bold');
ylabel('S/N (dB)','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
legend('Joint (source 1)','Joint (source 2)','Traditional (source 1)','Traditional (source 2)','Location','Best');
text(-3*40/23,14+(65-60)*10/120,'b)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','left');

subplot(2,2,4);
plot([1:5*5],mses,'r','linewidth',3);
xlabel('Iteration NO #','Fontsize',12,'fontweight','bold');
ylabel('MSE','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(-3*25/23,8+(65-60)*8/120,'c)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','left');
print(gcf,'-depsc','-r300','fig4.eps');
% 
% ylabel('Time (s)','Fontsize',12,'fontweight','bold');
% title('Results (joint iter #5)','Fontsize',12,'fontweight','bold');


