clc;clear;close all;
load test_real6.mat 

% yc_snr(d1,d1b)
% yc_snr(d1,D11)
% yc_snr(d1,D33)
% 
% yc_snr(d2,d2b)
% yc_snr(d2,D22)
% yc_snr(d2,D44)

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
subplot(5,1,1:2);yc_imagesc(comp1(1:1000,:),98,1,x,t(1:1000));
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

subplot(5,1,3:4);yc_imagesc(comp2(1:1000,:),98,1,x,t(1:1000));
text(-50,-0.2,'b)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','left');
xticks(xts);
set(gca,'xticklabel',{'30','60','90','120'});
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

