function prs
%绘制共中心点地震记录图，按震中距排列，地震事件记录图亦可，并将理论S波到时和P波到时连成线
%可批量绘制很多不同CDP记录图

% CDP fold
cf = 'C:\Users\Administrator\Documents\MATLAB\bishe\data\CDP_qre';  %有很多CDP文件夹，每个文件夹里有待画到一张图上的sac文件
% CDP_plot fold
cp = 'C:\Users\Administrator\Documents\MATLAB\bishe\plot\CDP_plot';
% reserved time截取时间段
tr = 2000;
% sampling time
dt = 1;
% velocity model
pfad='prem.nd';       %可换其他模型，如'ak135.nd', 'iasp91.nd'


model=mkreadnd(pfad);   %读取速度模型
foldlist = dir(cf);
[en,em] = size(foldlist);
for i = 3:en
    filelist = dir(strcat(ef,'\',foldlist(i).name));
    [fn,fm] = size(filelist);
    for j = 3:fn
        clear pdata;
        sac=strcat(ef,'\',foldlist(i).name,'\',filelist(j).name);
        [SeisData,HdrData] = readsac_Simons(sac,0);
        dist(j-2) = HdrData.GCARC;
        To(j-2) = HdrData.O;       % 地震起始时间；
        T0(j-2) = min(mkttime('P',HdrData.GCARC,HdrData.EVDP,model));   %计算P波到时
        T1(j-2) = min(mkttime('S',HdrData.GCARC,HdrData.EVDP,model));   %计算S波到时
        pdata(:,1) = SeisData(To(j-2)/dt+1:(To(j-2)+tr)/dt+1);   %截取数据从地震发生时刻开始：
        pdata = (pdata-SeisData(To(j-2)/dt+1))/(max(pdata)-SeisData(To(j-2)/dt+1))+dist;%调振幅使道间均衡
        partData(:,j-2)=pdata;
        t=0:dt:tr;
        h=plot(t,pdata,'k');
        xlabel('Time  /second');
        ylabel('Distance  /degree');
        set(gca,'YDir','reverse');
        hold on;
    end
    plot(T0,dist,'r',T1,dist,'b');  %将理论P到时（或S到时）连成线
    hold on;
    saveas(h,strcat(cp,'\',sprintf('%s',foldlist(i).name),'.fig'));   %保存记录图
    hlod off;
end

