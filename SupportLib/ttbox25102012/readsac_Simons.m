function [SeisData,HdrData]=readsac_Simons(filename,plotornot)
% [SeisData,HdrData]=READSAC(filename,plotornot)
% Reads in data saved from within SAC and plots or not (plotornot = 1: plot)
% SeisData is the seismic data
% HdrData contains information about the earthquake and receiver
% modified by FJS April 23th 1998
% last modified by Huajian Yao 2007

ppath=matlabpath;
    fid=fopen(filename,'rb');
    
if fid==-1
  error([ 'File ',filename,' does not exist in current path ',pwd]);
end
HdrFloats=fread(fid,70,'float32');
HdrNhdr=fread(fid,15,'int32');
HdrIhdr=fread(fid,20,'int32');
HdrLhdr=fread(fid,5,'int32');
HeaderStrings=str2mat(fread(fid,[8 24],'char'))';
SeisData=fread(fid,HdrNhdr(10),'float32');
fclose(fid);

HdrData=struct(...
  'NPTS',HdrNhdr(10),...                    % number of points
  'DELTA',HdrFloats(1),...                  % sampling time
  'SCALE',HdrFloats(4),...                  % amplitude scaling factor
  'B',HdrFloats(6),...                      % begin time of record
  'E',HdrFloats(7),...                      % end time of record
  'O',HdrFloats(8),...                      % event origin time (seconds relative to reference recording time)
  'T0',HdrFloats(11),...                    % timepick1(P)
  'T1',HdrFloats(12),...                    % timepick2(PP)
  'T2',HdrFloats(13),...                    % timepick3(PcP)
  'T3',HdrFloats(14),...                    % timepick4(S)
  'T4',HdrFloats(15),...                    % timepick5(SS)
  'T5',HdrFloats(16),...                    % timepick6(ScS)
  'KZYEAR',HdrNhdr(1),...                   % year
  'KZJDAY',HdrNhdr(2),...                   % julian day
  'KZHOUR',HdrNhdr(3),...                   % hour
  'KZMIN',HdrNhdr(4),...                    % minute
  'KZSEC',HdrNhdr(5),...                    % second
  'KZMSEC',HdrNhdr(6),...                   % milisecond
  'KSTNM',deblank(HeaderStrings(1,:)),...   % station name
  'KCMPNM',deblank(HeaderStrings(21,:)),... % recording component
  'KNETWK',deblank(HeaderStrings(22,:)),... % station network      
  'KINST',deblank(HeaderStrings(24,:)),...  % generic name of recording instrument  
  'STLA',HdrFloats(32),...                  % station latitude
  'STLO',HdrFloats(33),...                  % station longitude
  'STEL',HdrFloats(34),...                  % station elevation
  'EVLA',HdrFloats(36),...                  % event latitude
  'EVLO',HdrFloats(37),...                  % event longitude
  'EVDP',HdrFloats(39),...                  % event depth
  'DIST',HdrFloats(51),...                  % epicentral distance in km
  'AZ',HdrFloats(52),...                    % azimuth
  'BAZ',HdrFloats(53),...                   % back azimuth
  'GCARC',HdrFloats(54));                   % epicentral distance in degrees between source and receiver

if plotornot==1
  plot(linspace(HdrData.B,HdrData.E,HdrData.NPTS),SeisData);  
  title([filename]);
  xlabel([ 'Time (s)']);
end