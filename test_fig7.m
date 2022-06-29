%% Reproducible script for Figure 7 (one hour?)

close all; clc;clear;

%% Please change the directory path
% requiring the DRR package
% https://github.com/chenyk1990/MATdrr
addpath(genpath('~/MATdrr'));
addpath(genpath('../subroutines'));

%% please download data from https://drive.google.com/file/d/1ge0Mn_SB4LUsVgOBvATh0iISwGQahKh4/view?usp=sharing
load yc_fieldsr.mat
%% in this dataset
%there are two sources data3d(:,:,1:60) and data3d(:,:,61:120)
%each source contains 120 shots
%there are 60 receivers
%

d1=data3d(:,:,1);
d2=data3d(:,:,61);

h1=1:120;
h2=1:120;

dt=0.004;
t=[0:1500-1]*dt;
nt=1500;
nx=120;
%% apply dbt_dithering
randn('state',202122);
shift1=floor(0.1*randn(1,size(d1,2))/dt);   % shift of data1 to data2
shift2=-shift1;                             % shift of data2 to data1

d1shift=dbt_dither(d1,shift1);
d2shift=dbt_dither(d2,shift2);

figure;
subplot(1,2,1);imagesc(h1,t,d1shift);
subplot(1,2,2);imagesc(h2,t,d2shift);

%% blend
d1b=d1+d2shift;
d2b=d2+d1shift;
figure;
subplot(1,2,1);imagesc(h1,t,d1b);
subplot(1,2,2);imagesc(h2,t,d2b);

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
delta0(find(abs(delta-delta0)>8))=delta(find(abs(delta-delta0)>8))+round(16*(rand(size(find(abs(delta-delta0)>8)))-0.5));%exclude the outliers
del00=delta0;    %initial guess (fixed)

% figure;plot(delta);hold on;plot(delta0);legend('Ground truth','Initial','Location','Best');

del=delta;      %ground truth
d1b=d1+dbt_dither(d2,del);
d2b=d2+dbt_dither(d1,-del);

figure;
subplot(1,2,1);dbt_imagesc(d1b);
subplot(1,2,2);dbt_imagesc(d2b);

%% blending for all receivers
data3db=zeros(size(data3d));
for ir=1:60
    data3db(:,:,ir)=data3d(:,:,ir)+dbt_dither(data3d(:,:,ir+60),del);
    data3db(:,:,ir+60)=data3d(:,:,ir+60)+dbt_dither(data3d(:,:,ir),-del);
end


%% first receiver
mask1=ones(size(d1));
mask1=dbt_mutterv(mask1,1,42,4.1);

mask2=ones(size(d2));
mask2=dbt_mutterv(mask2,61,42,4.1);

bd=d1b;
% delta0=delta;
%% first receiver
niters=[2,4,6,8,10];
d1b=data3db(:,:,1);
d2b=data3db(:,:,61);
d1=data3d(:,:,1);
d2=data3d(:,:,61);
bd=d1b;
for iter_out=1:5
    D1=zeros(nt,nx);
    D2=zeros(nt,nx);
    for iter=1:niters(iter_out)
        fprintf('\n Iter %d \n',iter);
        D1T=dbt_dither(D1,-delta0);
        D2T=dbt_dither(D2,delta0);
        D1u = D1 + 0.5*(d1b-(D1+D2T));    % updated model
        D2u = D2 + 0.5*(d2b-(D1T+D2));    % updated model
        %         D1= seislet_denoise_2d(D1u,dip,'ps',3,0.1,2,3);
        %         D2= seislet_denoise_2d(D2u,dip2,'ps',3,0.1,2,3);
        D1u=D1u.*mask1;
        D2u=D2u.*mask2;
        D1=drr3d_win(D1u,0,80,0.004,2,4,0,100,20,1,0.5,0.5,0.5);
        D2=drr3d_win(D2u,0,80,0.004,2,4,0,100,20,1,0.5,0.5,0.5);
        D1=D1.*mask1;
        D2=D2.*mask2;
        
        D1=D1+(d1b-D1).*(1-mask2).*mask1;
        D2=D2+(d2b-D2).*(1-mask1).*mask2;
        
        snr2(iter)=10*log10(sum(sum(d1.*d1))/sum(sum((d1-D1).*(d1-D1))));
        snr22(iter)=10*log10(sum(sum(d2.*d2))/sum(sum((d2-D2).*(d2-D2))));
        figure(9);
        subplot(1,2,1);dbt_imagesc(D1);
        subplot(1,2,2);dbt_imagesc(D2);
        fprintf('iter2=%d,iter1=%d,SNR=%g,SNR2=%g\n',iter_out,iter,snr2(iter),snr22(iter));
    end
    
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
        left=dbt_bandpass(left,0.004,0,8,6,6,0,0);
        A=zeros(n1*n2,n2);
        for i2=1:n2
            A(1+(i2-1)*n1:i2*n1,i2)=left(:,i2);
        end
        bd0=D1+dbt_dither(D2,del0);
        %           bd0=d1+dbt_dither(d2,del0);
        right=bd0-bd;
        %         right=bd0(1:10:end,:)-bd(1:10:end,:);
        right=dbt_bandpass(right,0.004,0,8,6,6,0,0);
        dd=inv(A'*A)*(A'*right(:));
        %         randn('state',2000+iter); %random update does not work
        %         dd=randn(size(dd))*0.5;   %random update does not work
        del1=round(del0+1*dd);
    end
    delta0=del1;
    fprintf('iter2=%d,SNR=%g->%g,MSE=%g->%g\n',iter_out,dbt_snr(del(:),del00(:)),dbt_snr(del(:),del1(:)),sum((del(:)-del00(:)).^2)/nx,sum((del(:)-del1(:)).^2)/nx);
end

%% for all receivers
data3ddb=zeros(size(data3db));
niter=10;
for ir=1:60
    mask1=ones(size(d1));
    mask1=dbt_mutterv(mask1,ir,42,4.1);
    
    mask2=ones(size(d2));
    mask2=dbt_mutterv(mask2,60+ir,42,4.1);
    d1b=data3db(:,:,ir);
    d2b=data3db(:,:,60+ir);
    d1=data3d(:,:,ir);
    d2=data3d(:,:,60+ir);
    
    D1=zeros(nt,nx);
    D2=zeros(nt,nx);
    for iter=1:niter
        %         fprintf('\n Iter %d \n',iter);
        D1T=dbt_dither(D1,-delta0);
        D2T=dbt_dither(D2,delta0);
        D1u = D1 + 0.5*(d1b-(D1+D2T));    % updated model
        D2u = D2 + 0.5*(d2b-(D1T+D2));    % updated model
        %         D1= seislet_denoise_2d(D1u,dip,'ps',3,0.1,2,3);
        %         D2= seislet_denoise_2d(D2u,dip2,'ps',3,0.1,2,3);
        D1u=D1u.*mask1;
        D2u=D2u.*mask2;
        D1=drr3d_win(D1u,0,80,0.004,2,4,0,100,20,1,0.5,0.5,0.5);
        D2=drr3d_win(D2u,0,80,0.004,2,4,0,100,20,1,0.5,0.5,0.5);
        D1=D1.*mask1;
        D2=D2.*mask2;
        
        D1=D1+(d1b-D1).*(1-mask2).*mask1;
        D2=D2+(d2b-D2).*(1-mask1).*mask2;
        
        snr2(iter)=10*log10(sum(sum(d1.*d1))/sum(sum((d1-D1).*(d1-D1))));
        snr22(iter)=10*log10(sum(sum(d2.*d2))/sum(sum((d2-D2).*(d2-D2))));
%         figure(9);
%         subplot(1,2,1);dbt_imagesc(D1);
%         subplot(1,2,2);dbt_imagesc(D2); 
    end
    data3ddb(:,:,ir)=D1;
    data3ddb(:,:,ir+60)=D2;
    fprintf('ir=%d,snr2=%g,snr22=%g\n',ir,snr2(niter),snr22(niter));
end


data3ddb0=zeros(size(data3db));
niter=10;
for ir=1:60
    mask1=ones(size(d1));
    mask1=dbt_mutterv(mask1,ir,42,4.1);
    
    mask2=ones(size(d2));
    mask2=dbt_mutterv(mask2,60+ir,42,4.1);
    d1b=data3db(:,:,ir);
    d2b=data3db(:,:,60+ir);
    d1=data3d(:,:,ir);
    d2=data3d(:,:,60+ir);
    
    D1=zeros(nt,nx);
    D2=zeros(nt,nx);
    for iter=1:niter
        %         fprintf('\n Iter %d \n',iter);
        D1T=dbt_dither(D1,-del00(:));
        D2T=dbt_dither(D2,del00(:));
        D1u = D1 + 0.5*(d1b-(D1+D2T));    % updated model
        D2u = D2 + 0.5*(d2b-(D1T+D2));    % updated model
        %         D1= seislet_denoise_2d(D1u,dip,'ps',3,0.1,2,3);
        %         D2= seislet_denoise_2d(D2u,dip2,'ps',3,0.1,2,3);
        D1u=D1u.*mask1;
        D2u=D2u.*mask2;
        D1=drr3d_win(D1u,0,80,0.004,2,4,0,100,20,1,0.5,0.5,0.5);
        D2=drr3d_win(D2u,0,80,0.004,2,4,0,100,20,1,0.5,0.5,0.5);
        D1=D1.*mask1;
        D2=D2.*mask2;
        
        D1=D1+(d1b-D1).*(1-mask2).*mask1;
        D2=D2+(d2b-D2).*(1-mask1).*mask2;
        
        snr2(iter)=10*log10(sum(sum(d1.*d1))/sum(sum((d1-D1).*(d1-D1))));
        snr22(iter)=10*log10(sum(sum(d2.*d2))/sum(sum((d2-D2).*(d2-D2))));
%         figure(9);
%         subplot(1,2,1);dbt_imagesc(D1);
%         subplot(1,2,2);dbt_imagesc(D2);
        
    end
    data3ddb0(:,:,ir)=D1;
    data3ddb0(:,:,ir+60)=D2;
    fprintf('ir2=%d,snr2=%g,snr22=%g\n',ir,snr2(niter),snr22(niter));
end



figure;plot(1:n2,del,'k-','linewidth',2);hold on;plot(1:n2,del00,'b-','linewidth',2);
plot(1:n2,del1,'r--','linewidth',2);
legend('True','Initial','Recovered','Location','Best');


for ir=10:10
    figure(1)
    subplot(1,2,1);dbt_imagesc([data3d(:,:,ir),data3db(:,:,ir),data3ddb(:,:,ir),data3db(:,:,ir)-data3ddb(:,:,ir)]);
    subplot(1,2,2);dbt_imagesc([data3d(:,:,ir+60),data3db(:,:,ir+60),data3ddb(:,:,ir+60),data3db(:,:,ir+60)-data3ddb(:,:,ir+60)]);
end

for is=10:10
    figure(2)
    subplot(1,2,1);dbt_imagesc([squeeze(data3d(:,is,1:60)),squeeze(data3db(:,is,1:60)),squeeze(data3ddb(:,is,1:60)),squeeze(data3db(:,is,1:60))-squeeze(data3ddb(:,is,1:60))]);
    subplot(1,2,2);dbt_imagesc([squeeze(data3d(:,is,61:120)),squeeze(data3db(:,is,61:120)),squeeze(data3ddb(:,is,61:120)),squeeze(data3db(:,is,61:120))-squeeze(data3ddb(:,is,61:120))]);
end

for ir=10:10
    figure(1)
    subplot(1,2,1);dbt_imagesc([data3d(:,:,ir),data3db(:,:,ir),data3ddb0(:,:,ir),data3db(:,:,ir)-data3ddb0(:,:,ir)]);
    subplot(1,2,2);dbt_imagesc([data3d(:,:,ir+60),data3db(:,:,ir+60),data3ddb0(:,:,ir+60),data3db(:,:,ir+60)-data3ddb0(:,:,ir+60)]);
end

for is=10:10
    figure(2)
    subplot(1,2,1);dbt_imagesc([squeeze(data3d(:,is,1:60)),squeeze(data3db(:,is,1:60)),squeeze(data3ddb0(:,is,1:60)),squeeze(data3db(:,is,1:60))-squeeze(data3ddb0(:,is,1:60))]);
    subplot(1,2,2);dbt_imagesc([squeeze(data3d(:,is,61:120)),squeeze(data3db(:,is,61:120)),squeeze(data3ddb0(:,is,61:120)),squeeze(data3db(:,is,61:120))-squeeze(data3ddb0(:,is,61:120))]);
end

for ir=10:10
    figure(1)
    subplot(2,2,1);dbt_imagesc([data3d(:,:,ir),data3db(:,:,ir),data3ddb(:,:,ir),data3db(:,:,ir)-data3ddb(:,:,ir)]);
    subplot(2,2,2);dbt_imagesc([data3d(:,:,ir+60),data3db(:,:,ir+60),data3ddb(:,:,ir+60),data3db(:,:,ir+60)-data3ddb(:,:,ir+60)]);
    subplot(2,2,3);dbt_imagesc([data3d(:,:,ir),data3db(:,:,ir),data3ddb0(:,:,ir),data3db(:,:,ir)-data3ddb0(:,:,ir)]);
    subplot(2,2,4);dbt_imagesc([data3d(:,:,ir+60),data3db(:,:,ir+60),data3ddb0(:,:,ir+60),data3db(:,:,ir+60)-data3ddb0(:,:,ir+60)]);

end


for is=10:10
    figure(2)
    subplot(2,2,1);dbt_imagesc([squeeze(data3d(:,is,1:60)),squeeze(data3db(:,is,1:60)),squeeze(data3ddb(:,is,1:60)),squeeze(data3db(:,is,1:60))-squeeze(data3ddb(:,is,1:60))]);
    subplot(2,2,2);dbt_imagesc([squeeze(data3d(:,is,61:120)),squeeze(data3db(:,is,61:120)),squeeze(data3ddb(:,is,61:120)),squeeze(data3db(:,is,61:120))-squeeze(data3ddb(:,is,61:120))]);
    subplot(2,2,3);dbt_imagesc([squeeze(data3d(:,is,1:60)),squeeze(data3db(:,is,1:60)),squeeze(data3ddb0(:,is,1:60)),squeeze(data3db(:,is,1:60))-squeeze(data3ddb0(:,is,1:60))]);
    subplot(2,2,4);dbt_imagesc([squeeze(data3d(:,is,61:120)),squeeze(data3db(:,is,61:120)),squeeze(data3ddb0(:,is,61:120)),squeeze(data3db(:,is,61:120))-squeeze(data3ddb0(:,is,61:120))]);
end


%% Figure 7
figure('units','normalized','Position',[0.2 0.4 0.4, 1],'color','w');
subplot(3,2,1);dbt_mada3d(data3d(:,:,1:60));title('Unblended','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);text(-30,0,'a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
xlabel('Shot #','Fontsize',12,'fontweight','normal','Rotation',-10,'horizontalalignment','center','verticalalignment','middle');
ylabel('Receiver #','Fontsize',12,'fontweight','normal','Rotation',20,'horizontalalignment','center','verticalalignment','middle');
zlabel('Time (s)','Fontsize',12,'fontweight','normal');
set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');
zlim([0,5]);

% subplot(3,4,2);dbt_mada3d(data3d(:,:,61:120));title('Source 2','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);
subplot(3,2,2);dbt_mada3d(data3db(:,:,1:60));title('Blended','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);text(-30,0,'b)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
xlabel('Shot #','Fontsize',12,'fontweight','normal','Rotation',-10,'horizontalalignment','center','verticalalignment','middle');
ylabel('Receiver #','Fontsize',12,'fontweight','normal','Rotation',20,'horizontalalignment','center','verticalalignment','middle');
zlabel('Time (s)','Fontsize',12,'fontweight','normal');
set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');
% ax = gca;
% ax.XAxisLocation='origin';
zlim([0,5]);

subplot(3,2,3);dbt_mada3d(data3ddb0(:,:,1:60));title('Deblended (traditional)','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);text(-30,0,'c)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
xlabel('Shot #','Fontsize',12,'fontweight','normal','Rotation',-10,'horizontalalignment','center','verticalalignment','middle');
ylabel('Receiver #','Fontsize',12,'fontweight','normal','Rotation',20,'horizontalalignment','center','verticalalignment','middle');
zlabel('Time (s)','Fontsize',12,'fontweight','normal');
set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');
zlim([0,5]);

% subplot(3,4,2);dbt_mada3d(data3d(:,:,61:120));title('Source 2','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);
subplot(3,2,4);dbt_mada3d(data3db(:,:,1:60)-data3ddb0(:,:,1:60));title('Noise (traditional)','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);text(-30,0,'d)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
xlabel('Shot #','Fontsize',12,'fontweight','normal','Rotation',-10,'horizontalalignment','center','verticalalignment','middle');
ylabel('Receiver #','Fontsize',12,'fontweight','normal','Rotation',20,'horizontalalignment','center','verticalalignment','middle');
zlabel('Time (s)','Fontsize',12,'fontweight','normal');
set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');
% ax = gca;
% ax.XAxisLocation='origin';
zlim([0,5]);

subplot(3,2,5);dbt_mada3d(data3ddb(:,:,1:60));title('Deblended (joint)','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);text(-30,0,'e)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
xlabel('Shot #','Fontsize',12,'fontweight','normal','Rotation',-10,'horizontalalignment','center','verticalalignment','middle');
ylabel('Receiver #','Fontsize',12,'fontweight','normal','Rotation',20,'horizontalalignment','center','verticalalignment','middle');
zlabel('Time (s)','Fontsize',12,'fontweight','normal');
set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');
zlim([0,5]);

% subplot(3,4,2);dbt_mada3d(data3d(:,:,61:120));title('Source 2','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);
subplot(3,2,6);dbt_mada3d(data3db(:,:,1:60)-data3ddb(:,:,1:60));title('Noise (joint)','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);text(-30,0,'f)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
xlabel('Shot #','Fontsize',12,'fontweight','normal','Rotation',-10,'horizontalalignment','center','verticalalignment','middle');
ylabel('Receiver #','Fontsize',12,'fontweight','normal','Rotation',20,'horizontalalignment','center','verticalalignment','middle');
zlabel('Time (s)','Fontsize',12,'fontweight','normal');
set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');
% ax = gca;
% ax.XAxisLocation='origin';
zlim([0,5]);

% Create textarrow
annotation(gcf,'textarrow',[0.782986111111111 0.739583333333334],...
[0.880478087649402 0.849933598937583],'String',{'Blending noise'},'linewidth',2,'fontweight','bold','color','r');
 
% Create textarrow
annotation(gcf,'textarrow',[0.305555555555556 0.262152777777778],...
[0.565737051792828 0.535192563081009],'String',{'Remaining noise'},'linewidth',2,'fontweight','bold','color','r');
 
% Create textarrow
annotation(gcf,'textarrow',[0.312500000000001 0.269097222222223],...
[0.863213811420982 0.832669322709163],'linewidth',2,'fontweight','bold','color','r');
 
% Create textarrow
annotation(gcf,'textarrow',[0.305555555555556 0.262152777777778],...
[0.264276228419654 0.233731739707835],'linewidth',2,'fontweight','bold','color','r');
 
% Create arrow
annotation(gcf,'arrow',[0.784722222222222 0.795138888888889],...
[0.877822045152722 0.847277556440903],'linewidth',2,'color','r');
 


print(gcf,'-depsc','-r300','fig7.eps');







