%% Reproducible script for Figure 6

clc;clear;close all;

%% Please change the directory path
% requiring the DRR package
% https://github.com/chenyk1990/MATdrr
addpath(genpath('~/MATdrr'));
addpath(genpath('subroutines'));

%% please download data from https://drive.google.com/file/d/1ge0Mn_SB4LUsVgOBvATh0iISwGQahKh4/view?usp=sharing
load yc_fieldsr.mat
%% in this dataset
%there are two sources data3d(:,:,1:60) and data3d(:,:,61:120)
%each source contains 120 shots
%there are 60 receivers
%

d1=data3d(:,:,1);
d2=data3d(:,:,60);
figure;
subplot(1,2,1);dbt_imagesc(d1);
subplot(1,2,2);dbt_imagesc(d2);

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

% dip=dbt_dip2d_i(d1,5,30,2,0.01, 1, 0.000001,[20,5,1],1);
% dip2=dbt_dip2d_i(d2,5,30,2,0.01, 1, 0.000001,[20,5,1],1);

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
% delta0(find(abs(delta-delta0)>12))=delta(find(abs(delta-delta0)>12))+round(24*(rand(size(find(abs(delta-delta0)>12)))-0.5));%exclude the outliers
% delta0(find(abs(delta-delta0)>4))=delta(find(abs(delta-delta0)>4))+round(8*(rand(size(find(abs(delta-delta0)>4)))-0.5));%exclude the outliers

del00=delta0;    %initial guess (fixed)

% figure;plot(delta);hold on;plot(delta0);legend('Ground truth','Initial','Location','Best');

del=delta;      %ground truth
d1b=d1+dbt_dither(d2,del);
d2b=d2+dbt_dither(d1,-del);

figure;
subplot(1,2,1);dbt_imagesc(d1b);
subplot(1,2,2);dbt_imagesc(d2b);


%% mask1
dd=[fliplr(d1),d1];
d1m=dbt_mutter(dd,120,42,505);
% figure;dbt_imagesc([dd,d1m,dd-d1m]);
mask1=ones(size(dd));
mask1=dbt_mutter(mask1,120,42,505);
mask1=mask1(:,121:end);
%mask2
mask2=ones(size(d2));
mask2=dbt_mutter(mask2,60,40,310);

% dd=[d2,fliplr(d2)];
% d2m=dbt_mutter(dd,120,42,505);
% figure;dbt_imagesc([dd,d2m,dd-d2m]);
% mask2=ones(size(dd));
% mask2=dbt_mutter(mask2,120,42,505);
% mask2=mask2(:,1:120);
figure;subplot(1,2,1);imagesc(mask1);subplot(1,2,2);imagesc(mask2);


% dtmp=d1b.*(1-fliplr(mask1)).*mask1;
% figure;dbt_imagesc([d1,dtmp,d1-dtmp]);

%%
% d2m=dbt_mutter(d2,60,40,310);
% figure;dbt_imagesc([d2,d2m,d2-d2m]);
%
% d2m=dbt_mutter(d2b,60,10,310);
% figure;dbt_imagesc([d2b,d2m,d2b-d2m]);
%
% %% test
% figure;dbt_imagesc([d1,d1.*mask1,d1-d1.*mask1]);
% figure;dbt_imagesc([d1b,d1b.*mask1,d1b-d1b.*mask1]);
% figure;dbt_imagesc([d2,d2.*mask2,d2-d2.*mask2]);
% figure;dbt_imagesc([d2b,d2b.*mask2,d2b-d2b.*mask2]);

bd=d1b;
% delta0=delta;
niters=[2,4,6,8,10];
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

%% maximum db iter2=5,iter1=10,SNR=11.1169,SNR2=10.6301



figure(10);
subplot(1,2,1);dbt_imagesc([d1,d1b,D1,d1b-D1]);
subplot(1,2,2);dbt_imagesc([d2,d2b,D2,d2b-D2]);

figure(10);
D11=D1+(d1b-D1).*(1-mask2).*mask1;
D22=D2+(d2b-D2).*(1-mask1).*mask2;
subplot(1,2,1);dbt_imagesc([d1,d1b,D11,d1b-D11]);
subplot(1,2,2);dbt_imagesc([d2,d2b,D22,d2b-D22]);


%% traditional methods
D1=zeros(nt,nx);
D2=zeros(nt,nx);
for iter=1:10
    fprintf('\n Iter %d \n',iter);
    D1T=dbt_dither(D1,-del00);
    D2T=dbt_dither(D2,del00);
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
%     figure(9);
%     subplot(1,2,1);dbt_imagesc(D1);
%     subplot(1,2,2);dbt_imagesc(D2);
    fprintf('iter2=%d,iter1=%d,SNR=%g,SNR2=%g\n',iter_out,iter,snr2(iter),snr22(iter));
end

D33=D1+(d1b-D1).*(1-mask2).*mask1;
D44=D2+(d2b-D2).*(1-mask1).*mask2;


% dbt_snr(d1,d1b)
% dbt_snr(d1,D11)
% dbt_snr(d1,D33)
% 
% dbt_snr(d2,d2b)
% dbt_snr(d2,D22)
% dbt_snr(d2,D44)


%% Figure 6
ngap=10;
comp1=[d1,zeros(n1,ngap),d1b,zeros(n1,ngap),D11,zeros(n1,ngap),d1b-D11,zeros(n1,ngap),d1-D11,zeros(n1,ngap),D33,zeros(n1,ngap),d1b-D33,zeros(n1,ngap),d1-D33]; 
comp2=[d2,zeros(n1,ngap),d2b,zeros(n1,ngap),D22,zeros(n1,ngap),d2b-D22,zeros(n1,ngap),d2-D22,zeros(n1,ngap),D44,zeros(n1,ngap),d2b-D44,zeros(n1,ngap),d2-D44]; 
t=[0:1500-1]*0.004;
x=1:size(comp1,2);
xts1=[30,60,90];
xts2=xts1+ngap+nx;
xts3=xts1+ngap*2+nx*2;
xts4=xts1+ngap*3+nx*3;
xts5=xts1+ngap*4+nx*4;
xts6=xts1+ngap*5+nx*5;
xts7=xts1+ngap*6+nx*6;
xts8=xts1+ngap*7+nx*7;
xts=[xts1,xts2,xts3,xts4,xts5,xts6,xts7,xts8];
figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(5,1,1:2);dbt_imagesc(comp1(1:1000,:),98,1,x,t(1:1000));
text(-50,-0.2,'a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','left');
xticks(xts);
set(gca,'xticklabel',{'30','60','90'});
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(nx/2,0.16,'Unblended','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(nx+ngap+nx/2,0.16,'Blended','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*2+nx/2,0.16,'Deblended (new)','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*3+nx/2,0.16,'Noise (new)','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*4+nx/2,0.16,'Error (new)','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*5+nx/2,0.16,'Deblended (old)','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*6+nx/2,0.16,'Noise (old)','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*7+nx/2,0.16,'Error (old)','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');

text(nx/2,4.35,'Shot NO #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(nx+ngap+nx/2,4.35,'Shot NO #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*2+nx/2,4.35,'Shot NO #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*3+nx/2,4.35,'Shot NO #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*4+nx/2,4.35,'Shot NO #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*5+nx/2,4.35,'Shot NO #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*6+nx/2,4.35,'Shot NO #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*7+nx/2,4.35,'Shot NO #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');

subplot(5,1,3:4);dbt_imagesc(comp2(1:1000,:),98,1,x,t(1:1000));
text(-50,-0.2,'b)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','left');
xticks(xts);
set(gca,'xticklabel',{'30','60','90'});
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(nx/2,0.16,'Unblended','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(nx+ngap+nx/2,0.16,'Blended','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*2+nx/2,0.16,'Deblended (new)','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*3+nx/2,0.16,'Noise (new)','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*4+nx/2,0.16,'Error (new)','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*5+nx/2,0.16,'Deblended (old)','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*6+nx/2,0.16,'Noise (old)','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*7+nx/2,0.16,'Error (old)','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(nx/2,4.35,'Shot NO #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(nx+ngap+nx/2,4.35,'Shot NO #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*2+nx/2,4.35,'Shot NO #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*3+nx/2,4.35,'Shot NO #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*4+nx/2,4.35,'Shot NO #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*5+nx/2,4.35,'Shot NO #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*6+nx/2,4.35,'Shot NO #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*7+nx/2,4.35,'Shot NO #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');

subplot(5,1,5);
plot(1:n2,del,'k-','linewidth',2);hold on;plot(1:n2,del00,'b-','linewidth',2);
plot(1:n2,del1,'r--','linewidth',2);
text(-7,70,'c)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','left');
legend('True','Initial','Recovered','Location','Best');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
xlabel('Shot NO #','Fontsize',10,'fontweight','bold');
ylabel('Time shifts','Fontsize',10,'fontweight','bold');

%% adding arrows
% Create ellipse
annotation(gcf,'ellipse',...
    [0.619722222222222 0.755644090305445 0.0583333333333335 0.104913678618858],...
    'Color',[1 0 0],...
    'LineWidth',2);

% Create textarrow
annotation(gcf,'textarrow',[0.634722222222224 0.641666666666667],...
    [0.6917919236118 0.757024327330259],'Color',[1 0 0],...
    'String',{'Strong remaining noise'},...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold');

% Create ellipse
annotation(gcf,'ellipse',...
    [0.652666666666667 0.386314941196333 0.032333333333333 0.0903054448871182],...
    'Color',[1 0 0],...
    'LineWidth',2);

% Create textarrow
annotation(gcf,'textarrow',[0.643055555555556 0.655555555555556],...
    [0.338138925294889 0.381389252948886],'Color',[1 0 0],...
    'String',{'Strong remaining noise'},...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold');

% Create textarrow
annotation(gcf,'textarrow',[0.779166666666667 0.772222222222222],...
    [0.352555701179554 0.429826347732704],'Color',[1 0 0],...
    'String',{'Signal leakage'},...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold');

% Create textarrow
annotation(gcf,'textarrow',[0.766666666666667 0.791666666666667],...
    [0.672754434424818 0.752293577981651],'Color',[1 0 0],...
    'String',{'Signal leakage'},...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold');

% Create textarrow
annotation(gcf,'textarrow',[0.861111111111112 0.861111111111111],...
    [0.329373045171869 0.378768020969856],'Color',[1 0 0],...
    'String',{'Large error'},...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold');

% Create textarrow
annotation(gcf,'textarrow',[0.861111111111111 0.861111111111111],...
    [0.689792442288514 0.7391874180865],'Color',[1 0 0],...
    'String',{'Large error'},...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold');

print(gcf,'-depsc','-r300','fig6.eps');

