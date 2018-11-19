fold='F:\Project\RFPreAnalysis\InputData';
key='200*';
FoldDir=dir(fullfile(fold,key));
Num=length(FoldDir);
fid=fopen(fullfile(fold,'Eventlist.txt'),'w');
for loni=1:Num
    fprintf(fid,'%s\n',FoldDir(loni).name);
end
fclose(fid);