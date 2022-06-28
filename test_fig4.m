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


