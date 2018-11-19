function partprs
%
%
% event fold
ef= 'C:\Users\Administrator\Documents\MATLAB\bishe\data\event_pre1';
% CDP fold
cf = 'C:\Users\Administrator\Documents\MATLAB\bishe\data\CDP_qre';
% CDP_plot fold
cp = 'C:\Users\Administrator\Documents\MATLAB\bishe\plot\evtpart_plot';
% phasepick time fold
tf = 'C:\Users\Administrator\Documents\MATLAB\bishe\picktime\event_pt';
%partData (TO-50:T0+50)fold with normolization and plus dist  
pf = 'C:\Users\Administrator\Documents\MATLAB\bishe\partData\evt_pd';
%part seisdata flod without normolization
psf='C:\Users\Administrator\Documents\MATLAB\bishe\partsData\evt_psd';
% reserved time
tr = 2000;
% sampling time
dt = 1;
% velocity model
pfad='prem.nd';
%forword and backword last time
st1=50;   st2=50;

mkdir(cp);
mkdir(pf);
mkdir(psf);
model=mkreadnd(pfad);
foldlist = dir(ef);
[en,em] = size(foldlist);
for i = 3:en
    clear tp;
    clear partData;
    clear partsData;
    clear dist;
    filelist = dir(strcat(ef,'\',foldlist(i).name));
    [fn,fm] = size(filelist);
    for j = 3:fn
        clear pdata;
        sac=strcat(ef,'\',foldlist(i).name,'\',filelist(j).name);
        [SeisData,HdrData] = readsac_Simons(sac,0);
        tp(j-2,1)=min(mkttime('P',HdrData.GCARC,HdrData.EVDP,model));
        dist(j-2) = HdrData.GCARC;
        To(j-2) = HdrData.O;
        if (tr+HdrData.O)<=(HdrData.NPTS-1)*dt
            pdata(:,1) = SeisData(To(j-2)/dt+1:(To(j-2)+tr)/dt+1);
            padata = (pdata-SeisData(To(j-2)/dt+1))/(max(pdata)-SeisData(To(j-2)/dt+1))+dist(j-2);
            jqdata=pdata((tp(j-2,1)-st1)/dt+1:(tp(j-2,1)+st2)/dt+1,1);
            jdata=padata((tp(j-2,1)-st1)/dt+1:(tp(j-2,1)+st2)/dt+1,1);
            partsData(:,j-2)=jqdata;
            partData(:,j-2)=jdata; %partData is the data from (T0-st1) to (T0+st2)
            t=0:dt:st1+st2;
            h=plot(t,jdata,'k');
            set(gca,'YDir','reverse');
            hold on;
        end
    end
    
    pdsize=size(partData);
    pdf=strcat(pf,'\',sprintf('%s',foldlist(i).name),'.dat');
    fid=fopen(pdf,'wb');
    fseek(fid,0,-1);
    fwrite(fid,pdsize,'float32');
    fwrite(fid,partData,'float32');
    fclose(fid);
    
    psdsize=size(partsData);
    psdf=strcat(psf,'\',sprintf('%s',foldlist(i).name),'.dat');
    fid=fopen(psdf,'wb');
    fseek(fid,0,-1);
    fwrite(fid,psdsize,'float32');
    fwrite(fid,partsData,'float32');
    fclose(fid);
    
    xlabel('Time  /second','fontsize',15);
    ylabel('Distance  /degree','fontsize',15);
    d1=min(dist)-round((max(dist)-min(dist))/10);
    d2=max(dist)+round((max(dist)-min(dist))/10);
    set(gca,'Box','on','fontsize',13, 'YLim',[d1 d2]);
    hold on;
    T0(1:fn-2)=st1;
    plot(T0,dist,'r');
    text(max(T0),max(dist)+1,'P');
    hold on;
    saveas(h,strcat(cp,'\',sprintf('%s',foldlist(i).name),'.fig'));
    hold off;
end