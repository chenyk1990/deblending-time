%Joint deblending and shot scheduling inversion (already very close)
%%works well when maximum/minimum time shifts are 8 samples.
close all; clc;clear;
addpath(genpath('~/chenyk/matlibcyk'));
% delta0=delta0+round(32*(rand(1,n2)-0.5));


load fieldsrcyk.mat
%% in this dataset
%there are two sources data3d(:,:,1:60) and data3d(:,:,61:120)
%each source contains 120 shots
%there are 60 receivers
%

d1=data3d(:,:,1);
d2=data3d(:,:,61);

figure;
subplot(1,2,1);yc_imagesc(d1);colormap(seis);
subplot(1,2,2);yc_imagesc(d2);colormap(seis);
%
%
% d3=squeeze(data3d(:,1,:));
% d4=squeeze(data3d(:,60,:));
% figure;
% subplot(1,2,1);yc_imagesc([d1,d3,d1-d3]);colormap(seis);
% subplot(1,2,2);yc_imagesc(d2);colormap(seis);
%

h1=1:120;
h2=1:120;

dt=0.004;
t=[0:1500-1]*dt;
nt=1500;
nx=120;
%% apply seisdithering
randn('state',202122);
shift1=floor(0.1*randn(1,size(d1,2))/dt);   % shift of data1 to data2
shift2=-shift1;                             % shift of data2 to data1

d1shift=seisdither(d1,shift1);
d2shift=seisdither(d2,shift2);

figure;
subplot(1,2,1);imagesc(h1,t,d1shift);
subplot(1,2,2);imagesc(h2,t,d2shift);



%% blend
d1b=d1+d2shift;
d2b=d2+d1shift;
figure;
subplot(1,2,1);imagesc(h1,t,d1b);
subplot(1,2,2);imagesc(h2,t,d2b);

% dip=yc_dip2d_i(d1,5,30,2,0.01, 1, 0.000001,[20,5,1],1);
% dip2=yc_dip2d_i(d2,5,30,2,0.01, 1, 0.000001,[20,5,1],1);

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
delta0=floor(yc_meanf(delta,2,1,2));
% delta0(find(abs(delta-delta0)>16))=delta(find(abs(delta-delta0)>16))+round(32*(rand(size(find(abs(delta-delta0)>16)))-0.5));%exclude the outliers
delta0(find(abs(delta-delta0)>8))=delta(find(abs(delta-delta0)>8))+round(16*(rand(size(find(abs(delta-delta0)>8)))-0.5));%exclude the outliers
% delta0(find(abs(delta-delta0)>12))=delta(find(abs(delta-delta0)>12))+round(24*(rand(size(find(abs(delta-delta0)>12)))-0.5));%exclude the outliers
% delta0(find(abs(delta-delta0)>4))=delta(find(abs(delta-delta0)>4))+round(8*(rand(size(find(abs(delta-delta0)>4)))-0.5));%exclude the outliers

del00=delta0;    %initial guess (fixed)

% figure;plot(delta);hold on;plot(delta0);legend('Ground truth','Initial','Location','Best');

del=delta;      %ground truth
d1b=d1+seisdither(d2,del);
d2b=d2+seisdither(d1,-del);

figure;
subplot(1,2,1);yc_imagesc(d1b);colormap(seis);
subplot(1,2,2);yc_imagesc(d2b);colormap(seis);

%% blending for all receivers
data3db=zeros(size(data3d));
for ir=1:60
    data3db(:,:,ir)=data3d(:,:,ir)+seisdither(data3d(:,:,ir+60),del);
    data3db(:,:,ir+60)=data3d(:,:,ir+60)+seisdither(data3d(:,:,ir),-del);
end


% for ir=1:60
%     figure(1)
%     subplot(1,2,1);yc_imagesc([data3d(:,:,ir),data3db(:,:,ir)]);colormap(seis);
%     subplot(1,2,2);yc_imagesc([data3d(:,:,ir+60),data3db(:,:,ir+60)]);colormap(seis);
%     pause;
% end
%
% for is=1:120
%     figure(1)
%     subplot(1,2,1);yc_imagesc([squeeze(data3d(:,is,1:60)),squeeze(data3db(:,is,1:60))]);colormap(seis);
%     subplot(1,2,2);yc_imagesc([squeeze(data3d(:,is,61:120)),squeeze(data3db(:,is,61:120))]);colormap(seis);
%     pause;
% end

%% generate mask
% for ir=1:60
%     mask1=ones(size(d1));
%     mask1=yc_mutterv(mask1,ir,42,4.1);
%
%     mask2=ones(size(d2));
%     mask2=yc_mutterv(mask2,60+ir,42,4.1);
%
%     figure(1)
%     subplot(1,2,1);yc_imagesc([data3d(:,:,ir),data3db(:,:,ir),data3db(:,:,ir).*mask1,data3d(:,:,ir)-data3d(:,:,ir).*mask1]);colormap(seis);
%     subplot(1,2,2);yc_imagesc([data3d(:,:,ir+60),data3db(:,:,ir+60),data3db(:,:,ir+60).*mask2,data3d(:,:,ir+60)-data3d(:,:,ir+60).*mask2]);colormap(seis);
%     pause;
% end

load data3ddb.mat
load data3ddb0.mat

figure('units','normalized','Position',[0.2 0.4 0.4, 1],'color','w');
subplot(3,2,1);yc_mada3d(data3d(:,:,1:60));title('Unblended','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);text(-30,0,'a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
xlabel('Shot #','Fontsize',12,'fontweight','normal','Rotation',-10,'horizontalalignment','center','verticalalignment','middle');
ylabel('Receiver #','Fontsize',12,'fontweight','normal','Rotation',20,'horizontalalignment','center','verticalalignment','middle');
zlabel('Time (s)','Fontsize',12,'fontweight','normal');
set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');
zlim([0,5]);

% subplot(3,4,2);yc_mada3d(data3d(:,:,61:120));title('Source 2','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);
subplot(3,2,2);yc_mada3d(data3db(:,:,1:60));title('Blended','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);text(-30,0,'b)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
xlabel('Shot #','Fontsize',12,'fontweight','normal','Rotation',-10,'horizontalalignment','center','verticalalignment','middle');
ylabel('Receiver #','Fontsize',12,'fontweight','normal','Rotation',20,'horizontalalignment','center','verticalalignment','middle');
zlabel('Time (s)','Fontsize',12,'fontweight','normal');
set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');
% ax = gca;
% ax.XAxisLocation='origin';
zlim([0,5]);

subplot(3,2,3);yc_mada3d(data3ddb0(:,:,1:60));title('Deblended (traditional)','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);text(-30,0,'c)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
xlabel('Shot #','Fontsize',12,'fontweight','normal','Rotation',-10,'horizontalalignment','center','verticalalignment','middle');
ylabel('Receiver #','Fontsize',12,'fontweight','normal','Rotation',20,'horizontalalignment','center','verticalalignment','middle');
zlabel('Time (s)','Fontsize',12,'fontweight','normal');
set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');
zlim([0,5]);

% subplot(3,4,2);yc_mada3d(data3d(:,:,61:120));title('Source 2','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);
subplot(3,2,4);yc_mada3d(data3db(:,:,1:60)-data3ddb0(:,:,1:60));title('Noise (traditional)','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);text(-30,0,'d)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
xlabel('Shot #','Fontsize',12,'fontweight','normal','Rotation',-10,'horizontalalignment','center','verticalalignment','middle');
ylabel('Receiver #','Fontsize',12,'fontweight','normal','Rotation',20,'horizontalalignment','center','verticalalignment','middle');
zlabel('Time (s)','Fontsize',12,'fontweight','normal');
set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');
% ax = gca;
% ax.XAxisLocation='origin';
zlim([0,5]);

subplot(3,2,5);yc_mada3d(data3ddb(:,:,1:60));title('Deblended (joint)','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);text(-30,0,'e)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
xlabel('Shot #','Fontsize',12,'fontweight','normal','Rotation',-10,'horizontalalignment','center','verticalalignment','middle');
ylabel('Receiver #','Fontsize',12,'fontweight','normal','Rotation',20,'horizontalalignment','center','verticalalignment','middle');
zlabel('Time (s)','Fontsize',12,'fontweight','normal');
set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');
zlim([0,5]);

% subplot(3,4,2);yc_mada3d(data3d(:,:,61:120));title('Source 2','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);
subplot(3,2,6);yc_mada3d(data3db(:,:,1:60)-data3ddb(:,:,1:60));title('Noise (joint)','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);text(-30,0,'f)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
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

% subplot(3,4,4);yc_mada3d(data3db(:,:,61:120));title('Source 2','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);
% 
% subplot(3,4,5);yc_mada3d(data3ddb0(:,:,1:60));title('Source 1','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);text(-50,-8,'c)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% subplot(3,4,6);yc_mada3d(data3ddb0(:,:,61:120));title('Source 2','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);
% subplot(3,4,7);yc_mada3d(data3db(:,:,1:60)-data3ddb0(:,:,1:60));title('Source 1','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);text(-50,-8,'d)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% subplot(3,4,8);yc_mada3d(data3db(:,:,61:120)-data3ddb0(:,:,61:120));title('Source 2','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);
% 
% subplot(3,4,9);yc_mada3d(data3ddb(:,:,1:60));title('Source 1','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);text(-50,-8,'e)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% subplot(3,4,10);yc_mada3d(data3ddb(:,:,61:120));title('Source 2','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);
% subplot(3,4,11);yc_mada3d(data3db(:,:,1:60)-data3ddb(:,:,1:60));title('Source 1','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);text(-50,-8,'f)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% subplot(3,4,12);yc_mada3d(data3db(:,:,61:120)-data3ddb(:,:,61:120));title('Source 2','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);
% 







