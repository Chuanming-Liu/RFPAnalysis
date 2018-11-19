function prs
%���ƹ����ĵ�����¼ͼ�������о����У������¼���¼ͼ��ɣ���������S����ʱ��P����ʱ������
%���������ƺܶ಻ͬCDP��¼ͼ

% CDP fold
cf = 'C:\Users\Administrator\Documents\MATLAB\bishe\data\CDP_qre';  %�кܶ�CDP�ļ��У�ÿ���ļ������д�����һ��ͼ�ϵ�sac�ļ�
% CDP_plot fold
cp = 'C:\Users\Administrator\Documents\MATLAB\bishe\plot\CDP_plot';
% reserved time��ȡʱ���
tr = 2000;
% sampling time
dt = 1;
% velocity model
pfad='prem.nd';       %�ɻ�����ģ�ͣ���'ak135.nd', 'iasp91.nd'


model=mkreadnd(pfad);   %��ȡ�ٶ�ģ��
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
        To(j-2) = HdrData.O;       % ������ʼʱ�䣻
        T0(j-2) = min(mkttime('P',HdrData.GCARC,HdrData.EVDP,model));   %����P����ʱ
        T1(j-2) = min(mkttime('S',HdrData.GCARC,HdrData.EVDP,model));   %����S����ʱ
        pdata(:,1) = SeisData(To(j-2)/dt+1:(To(j-2)+tr)/dt+1);   %��ȡ���ݴӵ�����ʱ�̿�ʼ��
        pdata = (pdata-SeisData(To(j-2)/dt+1))/(max(pdata)-SeisData(To(j-2)/dt+1))+dist;%�����ʹ�������
        partData(:,j-2)=pdata;
        t=0:dt:tr;
        h=plot(t,pdata,'k');
        xlabel('Time  /second');
        ylabel('Distance  /degree');
        set(gca,'YDir','reverse');
        hold on;
    end
    plot(T0,dist,'r',T1,dist,'b');  %������P��ʱ����S��ʱ��������
    hold on;
    saveas(h,strcat(cp,'\',sprintf('%s',foldlist(i).name),'.fig'));   %�����¼ͼ
    hlod off;
end

