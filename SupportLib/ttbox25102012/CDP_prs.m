function CDP_prs
%plot
%write picktime T0(P) into the dat file
% event fold
ef= 'C:\Users\Administrator\Documents\MATLAB\bishe\data\CDP_pre_1';
% CDP fold
cf = 'C:\Users\Administrator\Documents\MATLAB\bishe\data\CDP_qre';
% CDP_plot fold
cp = 'C:\Users\Administrator\Documents\MATLAB\bishe\plot\CDP_plot';
% phasepick time fold
tf = 'C:\Users\Administrator\Documents\MATLAB\bishe\picktime\CDP_pt';
%O starting data flod with normolization and plus dist
odf ='C:\Users\Administrator\Documents\MATLAB\bishe\OstarData\CDP_od';
%O starting seisdata fold without normalization
osdf='C:\Users\Administrator\Documents\MATLAB\bishe\OseisData\CDP_osd';
%dist fold
distf ='C:\Users\Administrator\Documents\MATLAB\bishe\dist\CDP_dist';
% raylength fold
rayf='C:\Users\Administrator\Documents\MATLAB\bishe\raylength\CDP_rayl';
% reserved time
tr = 2000;
% sampling time
dt = 1;
% velocity model
pfad='prem.nd';
%phase;   planet radius
phase='P';   rp=6371.004;


mkdir(cp);
mkdir(tf);
mkdir(odf);
mkdir(osdf);
mkdir(distf);
mkdir(rayf);
model=mkreadnd(pfad);
foldlist = dir(ef);
[en,em] = size(foldlist);
for i = 3:en
    clear tp;
    clear OstarData;
    clear Oseisdata;
    clear dist;
    clear rayl;
    filelist = dir(strcat(ef,'\',foldlist(i).name));
    [fn,fm] = size(filelist);
    k=0;
    for j = 3:fn
        clear pdata;
        clear padata;
        sac=strcat(ef,'\',foldlist(i).name,'\',filelist(j).name);
        [SeisData,HdrData] = readsac_Simons(sac,0);
        if (tr+HdrData.O)<=(HdrData.NPTS-1)*dt
            [tt,p]=mkttime(phase,HdrData.GCARC,HdrData.EVDP,model);
            np=find(tt==min(tt));
            tp(j-2-k,1)=min(tt);
            [d,segx,segz,segtyp,resp]=mkx4p(phase,HdrData.EVDP,p(np),model);
            rayl(j-2-k,1)=mkraylength(segx,segz,rp);
            %Wtpicks(tp(j-2-k,1),tp(j-2-k,2),tp(j-2-k,3),tp(j-2-k,4),tp(j-2-k,5),tp(j-2-k,6),sac);
            dist(j-2-k) = HdrData.GCARC;
            To(j-2-k) = HdrData.O;
            pdata(:,1) = SeisData(To(j-2-k)/dt+1:(To(j-2-k)+tr)/dt+1);
            padata = (pdata-SeisData(To(j-2-k)/dt+1))/(max(pdata)-SeisData(To(j-2-k)/dt+1))+dist(j-2-k);
            OseisData(:,j-2-k)=pdata;
            OstarData(:,j-2-k)=padata;
            t=0:dt:tr;
            h=plot(t,padata,'k');
            set(gca,'YDir','reverse');
            hold on;
        else
            k=k+1;
        end
    end
    
    tpsize=size(tp);
    odsize=size(OstarData);
    dsize=size(dist);
    osdsize=size(OseisData);
    rlsize=size(rayl);
    
    tpf=strcat(tf,'\',sprintf('%s',foldlist(i).name),'.dat');
    fid=fopen(tpf,'wb');
    fseek(fid,0,-1);
    fwrite(fid,tpsize,'float32');
    fwrite(fid,tp,'float32');
    fclose(fid);
    
    osf=strcat(odf,'\',sprintf('%s',foldlist(i).name),'.dat');
    fid=fopen(osf,'wb');
    fseek(fid,0,-1);
    fwrite(fid,odsize,'float32');
    fwrite(fid,OstarData,'float32');
    fclose(fid);
    
    disf=strcat(distf,'\',sprintf('%s',foldlist(i).name),'.dat');
    fid=fopen(disf,'wb');
    fseek(fid,0,-1);
    fwrite(fid,dsize,'float32');
    fwrite(fid,dist,'float32');
    fclose(fid);
    
    osef=strcat(osdf,'\',sprintf('%s',foldlist(i).name),'.dat');
    fid=fopen(osef,'wb');
    fseek(fid,0,-1);
    fwrite(fid,osdsize,'float32');
    fwrite(fid,OseisData,'float32');
    fclose(fid);
    
    raylf=strcat(rayf,'\',sprintf('%s',foldlist(i).name),'.dat');
    fid=fopen(raylf,'wb');
    fseek(fid,0,-1);
    fwrite(fid,rlsize,'float32');
    fwrite(fid,rayl,'float32');
    fclose(fid);
    
    xlabel('Time  /second','fontsize',15);
    ylabel('Distance  /degree','fontsize',15);
    t1=t(1);  t2=t(tr)+50*dt;
    d1=min(dist)-round((max(dist)-min(dist))/10);
    d2=max(dist)+round((max(dist)-min(dist))/10);
    set(gca,'Box','on','fontsize',13, 'XLim', [t1 t2], 'YLim',[d1 d2]);
    hold on;
    T0=tp(:,1)';
    %T1=tp(:,2)';
    %T2=tp(:,3)';
    %T3=tp(:,4)';
    %T4=tp(:,5)';
    %T5=tp(:,6)';
    T01=T0-50;
    T02=T0+50;
    plot(T0,dist,'r',T01,dist,'b',T02,dist,'b');
    %T1,dist,'r',T2,dist,'r',T3,dist,'b',T4,dist,'b',T5,dist,'b',
    text(max(T0),max(dist)+1,'P');
    %text(max(T1),max(dist)+1,'PP');
    %text(max(T2),max(dist)+1,'PcP');
    %text(max(T3),max(dist)+1,'S');
    %text(max(T4),max(dist)+1,'SS');
    %text(max(T5),max(dist)+1,'ScS');
    hold on;
    saveas(h,strcat(cp,'\',sprintf('%s',foldlist(i).name),'.fig'));
    hold off;
end


function Wtpicks(T0,T1,T2,T3,T4,T5,filename)
% write picked phase'times into the SAC files
% T0(P);T1(PP);T2(PcP);T3(S);T4(SS);T5(ScS);

fid=fopen(filename,'rb+');

if fid==-1
    error([ 'File ',filename,' does not exist in current path ',pwd]);
end
HdrFloats=fread(fid,70,'float32');
HdrNhdr=fread(fid,15,'int32');
HdrIhdr=fread(fid,20,'int32');
HdrLhdr=fread(fid,5,'int32');
HeaderStrings=str2mat(fread(fid,[8 24],'char'))';
SeisData=fread(fid,HdrNhdr(10),'float32');

HdrFloats(11)=T0+HdrData.O;
HdrFloats(12)=T1+HdrData.O;
HdrFloats(13)=T2+HdrData.O;
HdrFloats(14)=T3+HdrData.O;
HdrFloats(15)=T4+HdrData.O;
HdrFloats(16)=T5+HdrData.O;
fseek(fid,0,-1);
fwrite(fid,HdrFloats,'float32');
fclose(fid);

