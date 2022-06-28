% script for fig3
close all; clc;clear;
addpath(genpath('~/chenyk/matlibcyk'));

% save test_syn7.mat 
% save test_syn7new.mat 
% load test_syn7.mat
load test_syn7new.mat

% figure('units','normalized','Position',[0.2 0.4 0.6, 1],'color','w');
% for ii=1:4
% subplot(4,6,1+(ii-1)*6);imagesc(data1s(:,:,ii));caxis([-3,3]);colormap(seis);text(-10,-15,'a)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% subplot(4,6,2+(ii-1)*6);imagesc(data2s(:,:,ii));caxis([-3,3]);colormap(seis);
% subplot(4,6,3+(ii-1)*6);imagesc(d1b-data1s(:,:,ii));caxis([-3,3]);colormap(seis);text(-10,-15,'b)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% subplot(4,6,4+(ii-1)*6);imagesc(d2b-data2s(:,:,ii));caxis([-3,3]);colormap(seis);
% subplot(4,6,5+(ii-1)*6);imagesc(d1-data1s(:,:,ii));caxis([-3,3]);colormap(seis);text(-10,-15,'c)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% subplot(4,6,6+(ii-1)*6);imagesc(d2-data2s(:,:,ii));caxis([-3,3]);colormap(seis);
% end
% % print(gcf,'-depsc','-r300','fig_syn_iters.eps');
% 
% figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
% subplot(2,6,1);imagesc(data1s(:,:,5));caxis([-3,3]);colormap(seis);text(-10,-15,'a)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% subplot(2,6,2);imagesc(data2s(:,:,5));caxis([-3,3]);colormap(seis);
% subplot(2,6,3);imagesc(d1b-data1s(:,:,5));caxis([-3,3]);colormap(seis);text(-10,-15,'b)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% subplot(2,6,4);imagesc(d2b-data2s(:,:,5));caxis([-3,3]);colormap(seis);
% subplot(2,6,5);imagesc(d1-data1s(:,:,5));caxis([-3,3]);colormap(seis);text(-10,-15,'c)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% subplot(2,6,6);imagesc(d2-data2s(:,:,5));caxis([-3,3]);colormap(seis);
% 
% subplot(2,6,7);imagesc(D1);caxis([-3,3]);colormap(seis);text(-10,-15,'d)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% subplot(2,6,8);imagesc(D2);caxis([-3,3]);colormap(seis);
% subplot(2,6,9);imagesc(d1b-D1);caxis([-3,3]);colormap(seis);text(-10,-15,'e)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% subplot(2,6,10);imagesc(d2b-D2);caxis([-3,3]);colormap(seis);
% subplot(2,6,11);imagesc(d1-D1);caxis([-3,3]);colormap(seis);text(-10,-15,'f)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% subplot(2,6,12);imagesc(d2-D2);caxis([-3,3]);colormap(seis);
% % print(gcf,'-depsc','-r300','fig_syn_comp.eps');
% 
% figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
% subplot(2,1,1);
% plot(1:n2,del,'k-','linewidth',2);hold on;plot(1:n2,del00,'b-','linewidth',2);
% plot(1:n2,del1,'r--','linewidth',2); 
% xlabel('Shot NO');
% ylabel('Time shifts');
% legend('True','Initial','Recovered','Location','Best');
% 
% subplot(2,2,3);
% plot([1:8*5],snrs,'r','linewidth',2);hold on;plot([1:8*5],snrs2,'b','linewidth',2);
% plot([1:8],snr2,'g--','linewidth',2);
% plot([1:8],snr22,'c--','linewidth',2);
% xlabel('Iteration NO');
% ylabel('SNR');
% legend('Joint (source 1)','Joint (source 2)','Traditional (source 1)','Traditional (source 2)','Location','Best');
% 
% 
% subplot(2,2,4);
% plot([1:5*5],mses,'k','linewidth',2);
% xlabel('Iteration NO');
% ylabel('MSE');
% % print(gcf,'-depsc','-r300','fig_curves.eps');



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
set(gca,'xticklabel',{num2str(round(n2z/2))});
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
set(gca,'xticklabel',{num2str(round(n2z/2))});
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
set(gca,'xticklabel',{num2str(round(n2z/2))});
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
set(gca,'xticklabel',{num2str(round(n2z/2))});
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

