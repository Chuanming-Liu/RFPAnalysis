function Wttime
% If the starting time of the record is not not zero,the calcuated phase time 
% should add the parameter O;
%
% event fold
ef= 'C:\Users\Administrator\Documents\MATLAB\bishe\data\event_pre2';
% velocity model
pfad='prem.nd';


model=mkreadnd(pfad);
foldlist = dir(ef);
[en,em] = size(foldlist);
for i = 3:en
    filelist = dir(strcat(ef,'\',foldlist(i).name));
    [fn,fm] = size(filelist);
    for j = 3:fn
        clear pdata;
        sac=strcat(ef,'\',foldlist(i).name,'\',filelist(j).name);
        [SeisData,HdrData] = readsac_Simons(sac,0);
        t0=min(mkttime('P',HdrData.GCARC,HdrData.EVDP,model))+HdrData.O;
        t1=min(mkttime('PP',HdrData.GCARC,HdrData.EVDP,model))+HdrData.O;
        t2=min(mkttime('PcP',HdrData.GCARC,HdrData.EVDP,model))+HdrData.O;
        t3=min(mkttime('S',HdrData.GCARC,HdrData.EVDP,model))+HdrData.O;
        t4=min(mkttime('SS',HdrData.GCARC,HdrData.EVDP,model))+HdrData.O;
        t5=min(mkttime('ScS',HdrData.GCARC,HdrData.EVDP,model))+HdrData.O;
        Wtpicks(t0,t1,t2,t3,t4,t5,sac);
    end
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

HdrFloats(11)=T0;
HdrFloats(12)=T1;
HdrFloats(13)=T2;
HdrFloats(14)=T3;
HdrFloats(15)=T4;
HdrFloats(16)=T5;
fseek(fid,0,-1);
fwrite(fid,HdrFloats,'float32');
fclose(fid);
