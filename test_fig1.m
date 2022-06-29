%% Reproducible script for Figure 1
%
close all; clc;clear;

addpath(genpath('./subroutines'));
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

%% apply dbt_cseisdithering
randn('state',202122);
shift1=floor(0.1*randn(1,size(d1,2))/dt);   % shift of data1 to data2
shift2=-shift1;                             % shift of data2 to data1

d1shift=dbt_dither(d1,shift1);
d2shift=dbt_dither(d2,shift2);

%% blend
d1b=d1+d2shift;
d2b=d2+d1shift;
% figure;
% subplot(1,2,1);imagesc(h1,t,d1b);
% subplot(1,2,2);imagesc(h2,t,d2b);

[nt,nx]=size(d1);
t=[0:nt-1]*0.004;
x=1:nx;
% plotting figures
figure('units','normalized','Position',[0.2 0.4 0.4, 1],'color','w');
subplot(2,2,1);imagesc(x,t,d1);caxis([-2.5,2.5]);colormap(dbt_cseis);ylabel('Time (s)','Fontsize',12,'fontweight','bold');xlabel('Shot NO #','Fontsize',12,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');title('Source 1','Fontsize',12,'fontweight','bold');text(-10,-0.15,'a)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
subplot(2,2,2);imagesc(x,t,d2);caxis([-2.5,2.5]);colormap(dbt_cseis);ylabel('Time (s)','Fontsize',12,'fontweight','bold');xlabel('Shot NO #','Fontsize',12,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');title('Source 2','Fontsize',12,'fontweight','bold');text(-10,-0.15,'b)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
subplot(2,2,3);imagesc(x,t,d1b);caxis([-2.5,2.5]);colormap(dbt_cseis);ylabel('Time (s)','Fontsize',12,'fontweight','bold');xlabel('Shot NO #','Fontsize',12,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');title('Blended source 1','Fontsize',12,'fontweight','bold');text(-10,-0.15,'c)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
subplot(2,2,4);imagesc(x,t,d2b);caxis([-2.5,2.5]);colormap(dbt_cseis);ylabel('Time (s)','Fontsize',12,'fontweight','bold');xlabel('Shot NO #','Fontsize',12,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');title('Blended source 2','Fontsize',12,'fontweight','bold');text(-10,-0.15,'d)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
print(gcf,'-depsc','-r300','fig1.eps');


    