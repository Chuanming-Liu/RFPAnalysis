function phase_stack_qR
% stack the seisdata at a certain phase and cut off the Surface wave
%**************************************************************************
% event fold
ef= 'D:\lj\data\CDP_pre3_2';
% phasepick time fold
tf = 'D:\lj\picktime\CDP_pt\CDP_pre3_2\ScS';
%distant fold
distf='D:\lj\dist\CDP_dist\CDP_pre3_2';
%O starting seisdata fold without normalization
osdf='D:\lj\OseisData\CDP_osd\CDP_pre3_2';
% partstack data fold
pskf='D:\lj\partstack\CDP_psk\CDP_pre3_2\ScS';
% stack plot fold
spf='D:\lj\plot\stackplot\CDP_pre3_2';
% reserved time
tr=2700;
% the gradient range of the Surface wave
k1=32;    k2=50;
% sampling time
dt = 1;
% velocity model
pfad='prem.nd';
%phase
phase='ScS';   

mkdir(pskf);
mkdir(spf);
model=mkreadnd(pfad);
foldlist = dir(ef);
[en,em] = size(foldlist);
for i = 3:en
    clear partData tp dist OseisData partdata;
    tpf=strcat(tf,'\',sprintf('%s',foldlist(i).name),'.dat');
    fid=fopen(tpf,'rb');
    nst=fread(fid,1,'float32');
    nx=fread(fid,1,'float32');
    tp(1:nst)=fread(fid,nst,'float32');
    fclose(fid);
    
    disf=strcat(distf,'\',sprintf('%s',foldlist(i).name),'.dat');
    fid=fopen(disf,'rb');
    nx=fread(fid,1,'float32');
    nst=fread(fid,1,'float32');
    dist(1:nst)=fread(fid,nst,'float32');
    fclose(fid);

    tr1=min(mkttime(phase,0.03,0,model));
    tr1=round(tr1);
    tr2=tr-tr1;
    
    osef = strcat(osdf,'\',sprintf('%s',foldlist(i).name),'.dat');
    fid = fopen(osef,'rb');
    nt=fread(fid,1,'float32');
    nst=fread(fid,1,'float32');
    for j=1:nst
        OseisData(:,j)=fread(fid,nt,'float32');
        fprintf('the dist of file %s is %f\n',disf,dist(j));
        if dist(j)>0
        ts1=round(k1*dist(j)-100);  ts2=round(k2*dist(j)-50);
        OseisData(ts1/dt+1:ts2/dt+1,j)=0;
        end
    end
    fclose(fid);
    
    k=0;
    for j=1:nst
        n1=(round(tp(j))-tr1)/dt+1; n2=(round(tp(j))+tr2)/dt+1;
        if n1<=0
            n1=-n1+1;
            partData(1:n1,j)=0;
            partData(n1+1:n2+n1,j)=OseisData(1:n2,j);
        else
        partData(:,j)=OseisData(n1:n2,j);
        end
        
        if dist(j)<=40
            partdata(:,j-k)=partData(:,j);
        else
            k=k+1;
        end
    end
    
            
    np=length(partdata(1,:));     
    partdata=partdata';
    stackdata=sum(partdata)/np;
    sdata(:,i-2)=stackdata';
    
end
%reorder the data to make the grids at correct number
for i=1:6
    data(:,(i-1)*8+1:(i-1)*8+7)=sdata(:,(i-1)*8+2:(i-1)*8+8);
    data(:,(i-1)*8+8)=sdata(:,(i-1)*8+1);
end
for j=1:en-2
    Data(:,en-1-j)=data(:,j);
end

pssize=size(Data);
psf=strcat(pskf,'\','grid_t3.dat');
fid=fopen(psf,'wb');
fseek(fid,0,-1);
fwrite(fid,pssize,'float32');
fwrite(fid,Data,'float32');
fclose(fid);

filename=strcat(spf,'\','stack_ScS3.fig');
wigb1(Data,filename);
    

function wigb1(a,filename,scal,x,z,amx)
%WIGB1: Plot seismic data using wiggles
%
%  WIGB1(a,scal,x,z,amx) 
%
%  IN    a: seismic data
%        scale: multiple data by scale
%        x: x-axis;
%        z: vertical axis (time or depth)
%	 x and z are vectors with offset and time.
%
%	 If only 'a' is enter, 'scal,x,z,amn,amx' are decided automatically; 
%	 otherwise, 'scal' is a scalar; 'x, z' are vectors for annotation in 
%	 offset and time, amx are the amplitude range.
%
% Author:
% 	Xingong Li, Dec. 1995
% Changes: 
%	Jun11,1997: add amx
% 	May16,1997: updated for v5 - add 'zeros line' to background color
% 	May17,1996: if scal ==0, plot without scaling
% 	Aug6, 1996: if max(tr)==0, plot a line 

if nargin == 0, nx=10;nz=10; a = rand(nz,nx)-0.5; end;

[nz,nx]=size(a);

trmx= max(abs(a));
if (nargin <= 5); amx=mean(trmx);  end;
if (nargin <= 3); x=[1:nx]; z=[1:nz]; end;
if (nargin <= 2); scal =1; end;

if nx <= 1; disp(' ERR:PlotWig: nx has to be more than 1');return;end;

 % take the average as dx
	dx1 = abs(x(2:nx)-x(1:nx-1));
 	dx = median(dx1);

 dz=z(2)-z(1);
 xmx=max(max(a)); xmn=min(min(a)); 

 if scal == 0; scal=1; end;
 a = a * dx /amx; 
 a = a * scal;

 fprintf(' PlotWig: data range [%f, %f], plotted max %f \n',xmn,xmx,amx);
 
% set display range 
x1=min(x)-2.0*dx; x2=max(x)+2.0*dx;
z1=min(z)-dz; z2=max(z)+dz;
 
set(gca,'NextPlot','add','Box','on', ...
  'XLim', [x1 x2], ...
  'YDir','reverse', ...
  'YLim',[z1 z2]);
 
	fillcolor = [0 0 0];
	linecolor = [0 0 0];
	linewidth = 0.1;
    
    xlabel('grid number','FontName','Arial','FontSize',20);
    ylabel('t(s)','FontName','Arial','FontSize',20);
    
	z=z'; 	% input as row vector
	zstart=z(1);
	zend  =z(nz);
    set(gca,'FontName','Arial','FontSize',20)
for i=1:nx,
   
  if trmx(i) ~= 0;    % skip the zero traces
	tr=a(:,i); 	% --- one scale for all section
  	s = sign(tr) ;
  	i1= find( s(1:nz-1) ~= s(2:nz) );	% zero crossing points
	npos = length(i1);


	%12/7/97 
	zadd = i1 + tr(i1) ./ (tr(i1) - tr(i1+1)); %locations with 0 amplitudes
	aadd = zeros(size(zadd));

	[zpos,vpos] = find(tr >0);
	[zz,iz] = sort([zpos; zadd]); 	% indices of zero point plus positives
	aa = [tr(zpos); aadd];
	aa = aa(iz);

	% be careful at the ends
		if tr(1)>0, 	a0=0; z0=1.00;
		else, 		a0=0; z0=zadd(1);
		end;
		if tr(nz)>0, 	a1=0; z1=nz; 
		else, 		a1=0; z1=max(zadd);
		end;
			
	zz = [z0; zz; z1; z0];
 	aa = [a0; aa; a1; a0];
		

	zzz = zstart + zz*dz -dz;

	h=patch( aa+x(i) , zzz,  fillcolor);

	line( 'Color',[1 1 1],'EraseMode','background',  ...
         'Xdata', x(i)+[0 0], 'Ydata',[zstart zend]); % remove zero line

%'LineWidth',linewidth, ...
%12/7/97 'Xdata', x(i)+[0 0], 'Ydata',[z0 z1]*dz);	% remove zero line

	line( 'Color',linecolor,'EraseMode','background',  ...
	 'LineWidth',linewidth, ...
	 'Xdata', tr+x(i), 'Ydata',z);	% negatives line

   else % zeros trace
	line( 'Color',linecolor,'EraseMode','background',  ...
	 'LineWidth',linewidth, ...
         'Xdata', [x(i) x(i)], 'Ydata',[zstart zend]);
   end;
end;
hold on;
saveas(h,filename);