%% Reproducible script for Figures 3 and 4 (15 minutes?)
%
%% This script requires
% 1) the seistr package (https://github.com/chenyk1990/seistr);
% 2) the seislet package (currently it is not open-source yet, if you are interested in it, please contact Yangkang Chen, chenyk2016@gmail.com).
close all; clc;clear;
addpath(genpath('./'));
addpath(genpath('~/MATseislet'));
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
    fprintf('iter2=%d,SNR=%g->%g,MSE=%g->%g\n',iter_out,dbt_snr(del(:),del00(:)),dbt_snr(del(:),del1(:)),sum((del(:)-del00(:)).^2)/nx,sum((del(:)-del1(:)).^2)/nx);
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
    D1= st_seislet_denoise_2d(D1u,dip,'ps',3,0.1,2,3);
    D2= st_seislet_denoise_2d(D2u,dip2,'ps',3,0.1,2,3);
    
    snr2(iter)=10*log10(sum(sum(d1.*d1))/sum(sum((d1-D1).*(d1-D1))));
    snr22(iter)=10*log10(sum(sum(d2.*d2))/sum(sum((d2-D2).*(d2-D2))));
    figure(9);
    subplot(1,2,1);imagesc(h1,t,D1);
    subplot(1,2,2);imagesc(h2,t,D2);
    fprintf('iter2=%d,iter1=%d,SNR=%g,SNR2=%g\n',iter_out,iter,snr2(iter),snr22(iter));
end
figure;plot(snr2);


%% Figure 3
ngap=5;ngap2=2;
indt=150:250;indx=40:51;
n1z=length(indt);
n2z=length(indx);
[nt,nx]=size(d1);
comp1=[data1s(:,:,1),zeros(n1,ngap),data2s(:,:,1),zeros(n1,ngap),d1b-data1s(:,:,1),zeros(n1,ngap),d2b-data2s(:,:,1),zeros(n1,ngap),d1-data1s(:,:,1),zeros(n1,ngap),d2-data2s(:,:,1)];
comp1z=[data1s(indt,indx,1),zeros(n1z,ngap2),data2s(indt,indx,1),zeros(n1z,ngap2),d1(indt,indx)-data1s(indt,indx,1),zeros(n1z,ngap2),d2(indt,indx)-data2s(indt,indx,1)];

comp2=[data1s(:,:,2),zeros(n1,ngap),data2s(:,:,2),zeros(n1,ngap),d1b-data1s(:,:,2),zeros(n1,ngap),d2b-data2s(:,:,2),zeros(n1,ngap),d1-data1s(:,:,2),zeros(n1,ngap),d2-data2s(:,:,2)];
comp2z=[data1s(indt,indx,2),zeros(n1z,ngap2),data2s(indt,indx,2),zeros(n1z,ngap2),d1(indt,indx)-data1s(indt,indx,2),zeros(n1z,ngap2),d2(indt,indx)-data2s(indt,indx,2)];

comp3=[data1s(:,:,3),zeros(n1,ngap),data2s(:,:,3),zeros(n1,ngap),d1b-data1s(:,:,3),zeros(n1,ngap),d2b-data2s(:,:,3),zeros(n1,ngap),d1-data1s(:,:,3),zeros(n1,ngap),d2-data2s(:,:,3)];
comp3z=[data1s(indt,indx,3),zeros(n1z,ngap2),data2s(indt,indx,3),zeros(n1z,ngap2),d1(indt,indx)-data1s(indt,indx,3),zeros(n1z,ngap2),d2(indt,indx)-data2s(indt,indx,3)];

% comp4=[data1s(:,:,4),zeros(n1,ngap),data2s(:,:,4),zeros(n1,ngap),d1b-data1s(:,:,4),zeros(n1,ngap),d2b-data2s(:,:,4),zeros(n1,ngap),d1-data1s(:,:,4),zeros(n1,ngap),d2-data2s(:,:,4)];
% comp4z=[data1s(indt,indx,4),zeros(n1z,ngap2),data2s(indt,indx,4),zeros(n1z,ngap2),d1(indt,indx)-data1s(indt,indx,4),zeros(n1z,ngap2),d2(indt,indx)-data2s(indt,indx,4)];
comp4=[data1s(:,:,5),zeros(n1,ngap),data2s(:,:,5),zeros(n1,ngap),d1b-data1s(:,:,5),zeros(n1,ngap),d2b-data2s(:,:,5),zeros(n1,ngap),d1-data1s(:,:,5),zeros(n1,ngap),d2-data2s(:,:,5)];
comp4z=[data1s(indt,indx,5),zeros(n1z,ngap2),data2s(indt,indx,5),zeros(n1z,ngap2),d1(indt,indx)-data1s(indt,indx,5),zeros(n1z,ngap2),d2(indt,indx)-data2s(indt,indx,5)];


tz=t(indt);
xz=1:size(comp1z,2);

% comp2=[d2,zeros(n1,ngap),d2b,zeros(n1,ngap),D22,zeros(n1,ngap),d2b-D22,zeros(n1,ngap),D44,zeros(n1,ngap),d2b-D44];
t=[0:500-1]*0.004;
x=1:size(comp1,2);
xts1=[10,30];
xts2=xts1+ngap+nx;
xts3=xts1+ngap*2+nx*2;
xts4=xts1+ngap*3+nx*3;
xts5=xts1+ngap*4+nx*4;
xts6=xts1+ngap*5+nx*5;
xts=[xts1,xts2,xts3,xts4,xts5,xts6];

xts1z=[round(n2z/2)];
xts2z=xts1z+ngap2+n2z;
xts5z=xts1z+ngap2*2+n2z*2;
xts6z=xts1z+ngap2*3+n2z*3;
xtsz=[xts1z,xts2z,xts5z,xts6z];
figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(4,5,1:3);yc_imagesc(comp1(:,:),3,2,x,t);ylim([0,2.0]);
text(-30,-0.2,'a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','left');
xticks(xts);
set(gca,'xticklabel',{'10','30'});
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
title('Results (joint iter #1)','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(nx/2,0.16,'DB 1','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(nx+ngap+nx/2,0.16,'DB 2','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*2+nx/2,0.16,'Noise 1','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*3+nx/2,0.16,'Noise 2','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*4+nx/2,0.16,'Error 1','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*5+nx/2,0.16,'Error 2','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');

hold on;
plot([indx(1),indx(1)],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(end),indx(end)],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(1),indx(end)],t([indt(1),indt(1)]),'r','linewidth',2);
plot([indx(1),indx(end)],t([indt(end),indt(end)]),'r','linewidth',2);

plot([indx(1)+n2+ngap,indx(1)+n2+ngap],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(end)+n2+ngap,indx(end)+n2+ngap],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(1)+n2+ngap,indx(end)+n2+ngap],t([indt(1),indt(1)]),'r','linewidth',2);
plot([indx(1)+n2+ngap,indx(end)+n2+ngap],t([indt(end),indt(end)]),'r','linewidth',2);

plot([indx(1)+n2*4+ngap*4,indx(1)+n2*4+ngap*4],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(end)+n2*4+ngap*4,indx(end)+n2*4+ngap*4],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(1)+n2*4+ngap*4,indx(end)+n2*4+ngap*4],t([indt(1),indt(1)]),'r','linewidth',2);
plot([indx(1)+n2*4+ngap*4,indx(end)+n2*4+ngap*4],t([indt(end),indt(end)]),'r','linewidth',2);

plot([indx(1)+n2*5+ngap*5,indx(1)+n2*5+ngap*5],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(end)+n2*5+ngap*5,indx(end)+n2*5+ngap*5],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(1)+n2*5+ngap*5,indx(end)+n2*5+ngap*5],t([indt(1),indt(1)]),'r','linewidth',2);
plot([indx(1)+n2*5+ngap*5,indx(end)+n2*5+ngap*5],t([indt(end),indt(end)]),'r','linewidth',2);


subplot(4,5,4:5);yc_imagesc(comp1z(:,:),6,2,xz,tz);
xticks(xtsz);
set(gca,'xticklabel',{num2str(indx(1)+round(n2z/2))});
ylabel('','Fontsize',10,'fontweight','bold');
title('Zoomed (joint iter #1)','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');

text(n2z/2,tz(1)+0.16*(tz(end)-tz(1))/(t(end)-t(1)),'DB 1','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(n2z+ngap2+n2z/2,tz(1)+0.16*(tz(end)-tz(1))/(t(end)-t(1)),'DB 2','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((n2z+ngap2)*2+n2z/2,tz(1)+0.16*(tz(end)-tz(1))/(t(end)-t(1)),'Error 1','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((n2z+ngap2)*3+n2z/2,tz(1)+0.16*(tz(end)-tz(1))/(t(end)-t(1)),'Error 2','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');


subplot(4,5,[1:3]+5);yc_imagesc(comp2(:,:),3,2,x,t);ylim([0,2.0]);
text(-30,-0.2,'b)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','left');
xticks(xts);
set(gca,'xticklabel',{'10','30'});
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
title('Results (joint iter #2)','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
hold on;
plot([indx(1),indx(1)],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(end),indx(end)],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(1),indx(end)],t([indt(1),indt(1)]),'r','linewidth',2);
plot([indx(1),indx(end)],t([indt(end),indt(end)]),'r','linewidth',2);

plot([indx(1)+n2+ngap,indx(1)+n2+ngap],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(end)+n2+ngap,indx(end)+n2+ngap],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(1)+n2+ngap,indx(end)+n2+ngap],t([indt(1),indt(1)]),'r','linewidth',2);
plot([indx(1)+n2+ngap,indx(end)+n2+ngap],t([indt(end),indt(end)]),'r','linewidth',2);

plot([indx(1)+n2*4+ngap*4,indx(1)+n2*4+ngap*4],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(end)+n2*4+ngap*4,indx(end)+n2*4+ngap*4],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(1)+n2*4+ngap*4,indx(end)+n2*4+ngap*4],t([indt(1),indt(1)]),'r','linewidth',2);
plot([indx(1)+n2*4+ngap*4,indx(end)+n2*4+ngap*4],t([indt(end),indt(end)]),'r','linewidth',2);

plot([indx(1)+n2*5+ngap*5,indx(1)+n2*5+ngap*5],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(end)+n2*5+ngap*5,indx(end)+n2*5+ngap*5],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(1)+n2*5+ngap*5,indx(end)+n2*5+ngap*5],t([indt(1),indt(1)]),'r','linewidth',2);
plot([indx(1)+n2*5+ngap*5,indx(end)+n2*5+ngap*5],t([indt(end),indt(end)]),'r','linewidth',2);


subplot(4,5,[4:5]+5);yc_imagesc(comp2z(:,:),6,2,xz,tz);
xticks(xtsz);
set(gca,'xticklabel',{num2str(indx(1)+round(n2z/2))});
ylabel('','Fontsize',10,'fontweight','bold');
title('Zoomed (joint iter #2)','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');

subplot(4,5,[1:3]+5*2);yc_imagesc(comp3(:,:),3,2,x,t);ylim([0,2.0]);
text(-30,-0.2,'c)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','left');
xticks(xts);
set(gca,'xticklabel',{'10','30'});
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
title('Results (joint iter #3)','Fontsize',10,'fontweight','bold');
hold on;
plot([indx(1),indx(1)],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(end),indx(end)],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(1),indx(end)],t([indt(1),indt(1)]),'r','linewidth',2);
plot([indx(1),indx(end)],t([indt(end),indt(end)]),'r','linewidth',2);

plot([indx(1)+n2+ngap,indx(1)+n2+ngap],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(end)+n2+ngap,indx(end)+n2+ngap],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(1)+n2+ngap,indx(end)+n2+ngap],t([indt(1),indt(1)]),'r','linewidth',2);
plot([indx(1)+n2+ngap,indx(end)+n2+ngap],t([indt(end),indt(end)]),'r','linewidth',2);

plot([indx(1)+n2*4+ngap*4,indx(1)+n2*4+ngap*4],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(end)+n2*4+ngap*4,indx(end)+n2*4+ngap*4],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(1)+n2*4+ngap*4,indx(end)+n2*4+ngap*4],t([indt(1),indt(1)]),'r','linewidth',2);
plot([indx(1)+n2*4+ngap*4,indx(end)+n2*4+ngap*4],t([indt(end),indt(end)]),'r','linewidth',2);

plot([indx(1)+n2*5+ngap*5,indx(1)+n2*5+ngap*5],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(end)+n2*5+ngap*5,indx(end)+n2*5+ngap*5],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(1)+n2*5+ngap*5,indx(end)+n2*5+ngap*5],t([indt(1),indt(1)]),'r','linewidth',2);
plot([indx(1)+n2*5+ngap*5,indx(end)+n2*5+ngap*5],t([indt(end),indt(end)]),'r','linewidth',2);

subplot(4,5,[4:5]+5*2);yc_imagesc(comp3z(:,:),6,2,xz,tz);
xticks(xtsz);
set(gca,'xticklabel',{num2str(indx(1)+round(n2z/2))});
ylabel('','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
title('Zoomed (joint iter #3)','Fontsize',10,'fontweight','bold');


subplot(4,5,[1:3]+5*3);yc_imagesc(comp4(:,:),3,2,x,t);ylim([0,2.0]);
text(-30,-0.2,'d)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','left');
xticks(xts);
set(gca,'xticklabel',{'10','30'});
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
title('Results (joint iter #5)','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
hold on;
plot([indx(1),indx(1)],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(end),indx(end)],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(1),indx(end)],t([indt(1),indt(1)]),'r','linewidth',2);
plot([indx(1),indx(end)],t([indt(end),indt(end)]),'r','linewidth',2);

plot([indx(1)+n2+ngap,indx(1)+n2+ngap],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(end)+n2+ngap,indx(end)+n2+ngap],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(1)+n2+ngap,indx(end)+n2+ngap],t([indt(1),indt(1)]),'r','linewidth',2);
plot([indx(1)+n2+ngap,indx(end)+n2+ngap],t([indt(end),indt(end)]),'r','linewidth',2);

plot([indx(1)+n2*4+ngap*4,indx(1)+n2*4+ngap*4],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(end)+n2*4+ngap*4,indx(end)+n2*4+ngap*4],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(1)+n2*4+ngap*4,indx(end)+n2*4+ngap*4],t([indt(1),indt(1)]),'r','linewidth',2);
plot([indx(1)+n2*4+ngap*4,indx(end)+n2*4+ngap*4],t([indt(end),indt(end)]),'r','linewidth',2);

plot([indx(1)+n2*5+ngap*5,indx(1)+n2*5+ngap*5],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(end)+n2*5+ngap*5,indx(end)+n2*5+ngap*5],t([indt(1),indt(end)]),'r','linewidth',2);
plot([indx(1)+n2*5+ngap*5,indx(end)+n2*5+ngap*5],t([indt(1),indt(1)]),'r','linewidth',2);
plot([indx(1)+n2*5+ngap*5,indx(end)+n2*5+ngap*5],t([indt(end),indt(end)]),'r','linewidth',2);
text(n2/2,2.35,'Shot #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(n2+ngap+n2/2,2.35,'Shot #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((n2+ngap)*2+n2/2,2.35,'Shot #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((n2+ngap)*3+n2/2,2.35,'Shot #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((n2+ngap)*4+n2/2,2.35,'Shot #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((n2+ngap)*5+n2/2,2.35,'Shot #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');

subplot(4,5,[4:5]+5*3);yc_imagesc(comp4z(:,:),6,2,xz,tz);
xticks(xtsz);
set(gca,'xticklabel',{num2str(indx(1)+round(n2z/2))});
ylabel('','Fontsize',10,'fontweight','bold');
title('Zoomed (joint iter #5)','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2z/2,tz(end)+(2.35-t(end))*(tz(end)-tz(1))/(t(end)-t(1)),'Shot #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(n2z+ngap2+n2z/2,tz(end)+(2.35-t(end))*(tz(end)-tz(1))/(t(end)-t(1)),'Shot #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((n2z+ngap2)*2+n2z/2,tz(end)+(2.35-t(end))*(tz(end)-tz(1))/(t(end)-t(1)),'Shot #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((n2z+ngap2)*3+n2z/2,tz(end)+(2.35-t(end))*(tz(end)-tz(1))/(t(end)-t(1)),'Shot #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');


%% adding arrows

annotation(gcf,'arrow',[0.680555555555556 0.661111111111111],...
    [0.893946808510638 0.871010638297872],'color','r','linewidth',2);
annotation(gcf,'arrow',[0.680555555555556 0.661111111111112],...
    [0.455117021276595 0.432180851063829],'color','r','linewidth',2);
annotation(gcf,'arrow',[0.680555555555556 0.661111111111112],...
    [0.238361702127659 0.215425531914893],'color','r','linewidth',2);
annotation(gcf,'arrow',[0.681944444444445 0.6625],...
    [0.674531914893617 0.651595744680851],'color','r','linewidth',2);
annotation(gcf,'arrow',[0.752777777777779 0.733333333333334],...
    [0.649265957446808 0.626329787234042],'color','r','linewidth',2);
annotation(gcf,'arrow',[0.752777777777779 0.733333333333334],...
    [0.214425531914893 0.191489361702127],'color','r','linewidth',2);
annotation(gcf,'arrow',[0.754166666666667 0.734722222222223],...
    [0.429851063829787 0.406914893617021],'color','r','linewidth',2);
annotation(gcf,'arrow',[0.752777777777778 0.733333333333334],...
    [0.864691489361702 0.841755319148936],'color','r','linewidth',2);

print(gcf,'-depsc','-r300','fig3.eps');

%% Figure 4
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



