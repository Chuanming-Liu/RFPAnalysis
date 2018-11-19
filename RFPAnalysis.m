%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%    PRAnalysis     %%%%%%%%  %%
% The GUI is made for pretreatment for boady wave
% By Chuanming Liu, Aug 2015, at USTC
% Contact: lcm@mail.ustc.edu.cn
% Version 1.1.1
% Function to be added: 
% 1.screenshot
% 2.pick checkbox
% Verison 1.1.2
% Added:
% 1.calculate the P arrivel time automaticlly by ttbox(but too slow)
% 2.plot station map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                       Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = RFPAnalysis(varargin)
% RFPANALYSIS MATLAB code for RFPAnalysis.fig
%      RFPANALYSIS, by itself, creates a new RFPANALYSIS or raises the existing
%      singleton*.
%
%      H = RFPANALYSIS returns the handle to a new RFPANALYSIS or the handle to
%      the existing singleton
%
%      RFPANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RFPANALYSIS.M with the given input arguments.
%
%      RFPANALYSIS('Property','Value',...) creates a new RFPANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RFPAnalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RFPAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RFPAnalysis

% Last Modified by GUIDE v2.5 02-Mar-2016 21:09:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RFPAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @RFPAnalysis_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before RFPAnalysis is made visible.
function RFPAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RFPAnalysis (see VARARGIN)

% Choose default command line output for RFPAnalysis
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
handles.Tlength=45;
guidata(hObject, handles);
DefaultPath;

% UIWAIT makes RFPAnalysis wait for user response (see UIRESUME)
% uiwait(handles.PAFigure);

% --- Outputs from this function are returned to the command line.
function varargout = RFPAnalysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
axes1= handles.axes1;
AxesPlotSetting(axes1);

% --- Executes during object creation, after setting all properties.
function PAFigure_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PAFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
DataStructCreate;
cla;
% --- Executes when PAFigure is resized.
function PAFigure_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to PAFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Creat the DataStruture
function DataStructCreate
global TSdata PAdataInfo

global StationInfo SourceInfo WaveformInfo RecordInfo CrossCorrWaveInfo FilterInfo
PAdataInfo = struct('fname','',...
        'staname','',...
        'eventfile','',...
	    'handle',[],...
	    'quality',[],...
        'OrigWave',[],...
	    'wave',[],...
        'Fwave',[],...
        'GCDkm',[],...
        'GCDdeg',[],...
        'StartNum',[],...
        'EndNum',[],...
        'SENum',[],...
        'SampleT',[],...
        'SampleF',[],...
        'TimeL',[],...
        'PT',[],...
        'Lat',[],...
        'Lon',[]);
    
TSdata= struct('name','',...
	    'handle',[],...
	    'quality',[],...
	    'data',[]);
   
StationInfo = struct('Lon',0,...
	'Lat',0,...
    'SrcDepth',0,...
	'GCDkm',0,...
	'GCDdeg',0,...
	'Azim',0,...
    'SampleT',0,...
	'SampleF',0,...
    'Name','',...
    'Net','',...
    'Filename','');

RecordInfo = struct('YY',0,...
	'DD',0,...
	'HH',0,...
	'MM',0,...
	'SS',0,...
	'MS',0,...
	'DiffT',0,...
    'SampleT',0,...
	'SampleF',0,...
    'NumPt',0,...
    'Time',0,...
    'BT',[],...
    'ET',[],...
    'OT',[],...
    'PT',[]);
    
WaveformInfo = struct('DatZ',zeros(1,100),...
    'AmpZ',0);

SourceInfo = struct('Lon',0, ...
    'Lat',0, ...
    'Depth',0,...
    'YY',0,...
    'Month',0,...
    'Day',0,...
    'DD',0,...
	'HH',0,...
	'MM',0,...
	'SS',0,...
	'MS',0,...
    'Mag',0);

FilterInfo = struct('Mode',0,...
    'Domain',0,...
    'Window',0,...
    'CtrT',0,...
    'LowT',0,...
    'HighT',0,...
    'CtrF',0,...
    'LowF',0,...
    'HighF',0,...
    'SampleF',0,...
    'SampleT',0,...
    'Length',0,...
    'BandWidth',0,...
    'KaiserPara',1,...
    'GaussAlfa',2.5);

% Period range and interval for narrow band pass filtered cross correlation
CrossCorrWaveInfo = struct('StartT',0,...
    'EndT',0,...
    'StartF',0,...
    'EndF',0,...
    'DeltaT',0,...
    'DeltaF',0,...
    'NumCtrT',0,...
    'StartNum',0,...
    'EndNum',0,...
    'SENum',0,...
    'PointNum',0,...
    'WaveType',0,...
    'WinCode',0,...
    'HalfBand',0,...
    'wave1',zeros(1,100),...
    'wave2',zeros(1,100),...
    'group1',zeros(100,1),...
    'group2',zeros(100,1),...
    'VMin',0,...
    'VMax',0,...
    'DeltaV',0,...
    'PhasVImg',zeros(100,1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Main Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in StartingProcessing.
function StartingProcessing_Callback(hObject, eventdata, handles)
% hObject    handle to StartingProcessing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% the Main function
% the main cycly is station pair
% Initialization
% quanlity 1:good ;2:bad 0:vargin(event)
global Signal % 1:go on;2: event skip;3: exit
global eventlist EventlistDirectory  EvtsFold
global PAdata Increase
% Tip:eventlist.quanlity 0:vergin;1:good event; 2:bad event 
% Check
Begin=1;

if EventlistDirectory == 0
   Begin=0;
   eventlist_button_Callback(hObject, eventdata, handles);
   StartingProcessing_Callback(hObject, eventdata, handles);
else
   if isempty(eventlist)
      set(handles.MsgEdit,'String','Please Input the Eventlist File!');
      Begin=0;
   else
      event_Num=size(eventlist,2);
      if event_Num>0
         disp(['----Begin! EventNum=',num2str(event_Num),'----']);
         set(handles.EventNumAll,'String',num2str(event_Num));
      else
         Begin=0;
      end
   end
end

if Begin==0
   disp('Please check your input and output Path;then Push Start Button!')
else
% Initialization
   handles.index=[];
   guidata(hObject,handles); 
   handles.fig3=[];
   guidata(hObject,handles); 
   handles.increase=1;
   guidata(hObject,handles);
% error log
   disp(['The Outpat fold is:',EvtsFold]);
   errorlog=fullfile(EvtsFold,'SPAnalysis_log.txt');
   fid_log=fopen(errorlog,'a+');
   fprintf(fid_log,'\n **********************Start the PreAnalsysis***********************');
   fprintf(fid_log,'\n There are %s Events.',num2str(event_Num));
   fclose(fid_log);
   tic;
   index=0;
   Increase=1;
   Signal=1;
   while ((index+Increase)<=event_Num) && (Signal~=3)
       Signal=1;
       index=index+Increase;
       if index==0
           index=1;
       end
       clear PAdata
       handles.index=index;
       guidata(hObject,handles);
       handles=guidata(gcf);
       if ~isempty(handles.fig3)
           close(handles.fig3);
       end
       % judge Show all of the event or only show the unhandled and good
       % event 
       if get(handles.ShowAll,'Value')==1
          if eventlist(index).quality==0 || eventlist(index).quality==1 || eventlist(index).quality==2
              Begin=1;
          else
              Begin=0;
          end
       else
           if eventlist(index).quality==0 || eventlist(index).quality==1
              Begin=1;
           else 
               Begin=0;
           end
       end                        
       if Begin==1  
          set(handles.EventNumIndex,'String',num2str(index));
          word=['New earthquake event date: ',eventlist(index).foldname];
          set(handles.MsgEdit,'String',word);        
          event_fold_path=fullfile(EventlistDirectory,eventlist(index).foldname);
          if exist(event_fold_path,'dir')==0
             fid_log=fopen(errorlog,'a+');
             fprintf(fid_log,'\n the NO.%f event %s is not exist, go to the next event!',index,event_fold_path);
             disp(['the NO.',index,'event:',event_fold_path,'is not exist! Checke your evnetlist!']);
             fclose(fid_log);
             eventlist(index).quality=2;
             continue
          end
          % load the event filelist and buide src information 
          % the cycle to load every station record
          [sta,wave,record,source]=StationRecord_load(handles,event_fold_path,eventlist(index).foldname,errorlog);          
          if Signal==2 % mean null fold
             continue
          end  
          % PAdata sturcture establishment      
          PAdataSturcture(handles,sta,record,wave,eventlist(index).foldname,errorlog);
          if Signal==2 % mean all files do not contain effective data Or Ptime is negative
             continue
          end  
          % buide PAdata.wave 
          [ProcessIndex]=PAdataWaveBuide(hObject,eventdata,handles);      
          % load the saved file list if it exists.
          LoadFileList(handles,eventlist(index).foldname,errorlog);
          % Updata the information box
          UpdateInform(sta,source,eventlist(index).foldname,handles);
          % Filter Processing
          FilterProcess(handles);
          % choose to plot the orginal wave or filered wave
          if get(handles.FilterAuto,'Value')==1
             PlotFilterAxe1(handles);
          else
             % Plot the waveform on trace; Main figure
             PlotwaveformAxe1(handles);
          end          
          % Into CursorMode for the secondary control
          CursorMode(gcf);          
          uiwait(); 
          % Mk new event file or remove event 
          if Signal==3
             ExitCursorMode;
             break           
          end
       else
          word=['The NO.',num2str(index),'Bad event has been handled! Fold:',eventlist(index).foldname];
          disp(word);
          fid_log=fopen(errorlog,'a+');
          fprintf(fid_log,'\n %s',word);
          fclose(fid_log);
       end
   end
   set(handles.MsgEdit,'String','Data Processing Finished !');
   toc
end

%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%
%                    Loading the StationRecord files
function [sta,wave,record,source]=StationRecord_load(handles,event_fold_path,FoldName,error_log)
% read every station record of this event
% creat sta,wave,rcd_Z
% get the sta file list
global Signal eventlist
keyword=get(handles.CompKeyWord,'String');
if isempty(keyword)==1
   keyword='*BHZ*.SAC';
end
staDir=dir(fullfile(event_fold_path,keyword));
if isempty(staDir)
   Signal=2;
   eventlist(handles.index).quality=2;
   word=['The Event Fold Maybe Null! NO. ',num2str(handles.index),' Fold:',FoldName];
   disp(word);
   set(handles.MsgEdit,'String',word);
   fid_log=fopen(error_log,'a+');
   fprintf(fid_log,'\n----The Event fold maybe Null! NO.:%f \n %s',handles.index,FoldName);
   fclose(fid_log);
   sta=[];wave=[];record=[];source=[];
else
   fileNum=size(staDir,1);
   count=0;
   for loni=1:fileNum
       seisfile = fullfile(event_fold_path, staDir(loni).name);
       file_sta=dir(seisfile);
       if isempty(file_sta) 
          fid_log=fopen(error_log,'a+');
          fprintf(fid_log,'\n The station files is Null,\n %s',seisfile1);
          fclose(fid_log);
       elseif (file_sta.bytes< 3000) 
          fid_log=fopen(error_log,'a+');
          fprintf(fid_log,'\n The station files is too short:\n %s bytes is %f ',seisfile,file_sta.bytes);
          fclose(fid_log);
       else
          % load
          if Signal==1
             count=count+1;
             [sta(count),wave(count),record(count),source]=RdStaData(seisfile);
              sta(count).Filename=staDir(loni).name;
          else 
             break
          end
       end    
    end
end

% subfuction Rd one files
function [sta,wave,record,source]=RdStaData(seisfile)
global StationInfo WaveformInfo RecordInfo SourceInfo

sta = struct(StationInfo);
record = struct(RecordInfo);
wave = struct(WaveformInfo);
source = struct(SourceInfo);

S = readsac(seisfile);%readsac,readsac.simos
[sta, record, wave.DatZ, wave.AmpZ, source] = DataStructTrans(S,seisfile);


if (sta.GCDkm == -12345) || (sta.Azim == -12345) || isnan(sta.GCDkm ) || isnan(sta.Azim)
    sta.GCDkm = deg2km(distance([source.Lat,source.Lon],[sta.Lat,sta.Lon]));
    sta.Azim = azimuth([source.Lat,source.Lon],[sta.Lat,sta.Lon]);
end

% subfucntion transfrom sac file
function [station, record, ReSampleWave, ampcoef, source] = DataStructTrans(HdrData,seisfile)
global StationInfo SourceInfo RecordInfo
% HdrData contains information about the earthquake and receiver
% O:Earthquake Event origin time relative to reference absolut time
% B:Begining of the record time relative to reference absolut time
% P: manual pick P arrival time relative to reference absolut time
% record.DiffT = HdrData.O; Earthquake time- ref time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
station = struct(StationInfo);
source = struct(SourceInfo);
record = struct(RecordInfo);
station.Lat = HdrData.STLA;
station.Lon = HdrData.STLO;
station.GCDkm = HdrData.DIST;
station.GCDdeg = HdrData.GCARC;
station.Azim = HdrData.AZ;
station.SampleT = HdrData.DELTA;
station.SampleF = 1/station.SampleT;
station.Name = HdrData.KSTNM;
station.Net = HdrData.KNETWK;
station.SrcDepth = HdrData.EVDP;


record.YY = HdrData.NZYEAR;
record.DD = HdrData.NZJDAY;
record.HH = HdrData.NZHOUR;
record.MM = HdrData.NZMIN;
record.SS = HdrData.NZSEC;
record.MS = HdrData.NZMSEC;
record.DiffT = HdrData.O; 
record.NumPt = HdrData.NPTS;
record.SampleT = station.SampleT;
record.SampleF = station.SampleF;
record.OT = HdrData.O;
record.BT = HdrData.B;
record.ET = HdrData.E;
record.PT = HdrData.T0;
source.Depth = HdrData.EVDP;
source.Lat = HdrData.EVLA;
source.Lon = HdrData.EVLO;
source.Mag = HdrData.MAG;
% judge the  HdrData.B ==0 or not
% if abs(HdrData.B)>1
%    disp('The record time is not the reference time for this sac file!');
%    disp(seisfile);
% %  record.SS=record.SS+HdrData.B;
% %  error('The record time is not the reference time for this sac file!');
% end

if isnan(record.PT)
   Ptime=CalPTime(station.GCDdeg,source.Depth);
   record.PT=Ptime+record.OT;
end
record.PT=record.PT-record.BT;
% offset = mean(SeisData(1:record.NumPt));

SeisData(1:record.NumPt) = detrend(HdrData.DATA1(1:record.NumPt));
record.Time = (record.NumPt - 1)*station.SampleT;

% Resample the wave data
refSampleF=10;

if station.SampleF > refSampleF
    DecimateR = round(station.SampleF/refSampleF); 
    nn = floor(record.NumPt/DecimateR);
    if (nn*DecimateR+1) <= record.NumPt
        ReSampleWave = decimate(SeisData(1:nn*DecimateR+1)', DecimateR);
    else
        ReSampleWave = decimate([SeisData(1:nn*DecimateR)'; SeisData(nn*DecimateR)], DecimateR);
    end    

    % ReSampleWave = decimate(SeisData', DecimateR);
    station.SampleF = refSampleF;
    station.SampleT = 1/refSampleF;
    record.NumPt = length(ReSampleWave);
    record.SampleT = 1/refSampleF;
    record.SampleF = refSampleF;
else
    ReSampleWave = SeisData';
end
% ReSampleWave = SeisData';
ReSampleWave = ReSampleWave/max(ReSampleWave);
ampcoef = HdrData.SCALE; % amplitude scaling factor
record.Time = (record.NumPt - 1)*station.SampleT;

% calculate the P arrivel time by ttbox
function Ptime=CalPTime(GCARC,EVDP)
%  phase: string containing seismic phase name like 'P', 'S', 'ScS', 'PKPdf', etc.\
%          Phase names are case sensitive!
%  delta: epicentral distance [deg]
%  h: focal depth [km]
%  model: A structure describing the velocity distribution.
% velocity model
% Too slowly
pfad='prem.nd';
model=mkreadnd(pfad,'slient');
Ptime=mkttime('P',GCARC,EVDP,model);
%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%
% PAdata sturcture establishment
function PAdataSturcture(handles,sta,record,wave,eventfilename,errorlog)
% Tip PAdata.quality=1 good default ; PAdata.quality=2 bad
% Tip:eventlist.quanlity 0:vergin;1:good event; 2:bad event 
global PAdataInfo
global PAdata Signal eventlist
PAdata=struct(PAdataInfo);
StaNum=size(sta,2);
index=handles.index;
NullCount=0;
Effect=1;
for loni=1:StaNum
    if record(loni).PT<0 
       eventlist(handles.index).quality=2;
       disp(['P time is negative! NO.',num2str(handles.index),' Fold:',eventfilename,sta(loni).Filename]);
       fid_log=fopen(errorlog,'a+');
       fprintf(fid_log,'\n Error the NO.%f event %s,%s P time is negative!',index,eventfilename,sta(loni).Filename);
       fclose(fid_log);     
       Signal=2;
       break;
    end
    if isempty(wave(loni).DatZ) || all(isnan(wave(loni).DatZ))
       NullCount=NullCount+1;
    else   
       PAdata(Effect).eventfile=eventfilename;
       PAdata(Effect).fname=sta(loni).Filename;
       PAdata(Effect).staname=sta(loni).Name;
       PAdata(Effect).quality=1; % 1:good 2: bad
       PAdata(Effect).GCDkm=sta(loni).GCDkm;
       PAdata(Effect).GCDdeg=sta(loni).GCDdeg;
       PAdata(Effect).SampleT=record(loni).SampleT;
       PAdata(Effect).SampleF=record(loni).SampleF;
       PAdata(Effect).TimeL=record(loni).Time;
       PAdata(Effect).PT=record(loni).PT;
       PAdata(Effect).OrigWave=wave(loni).DatZ;
       PAdata(Effect).Lat=sta(loni).Lat;
       PAdata(Effect).Lon=sta(loni).Lon;
       Effect=Effect+1;
    end
end
% if all files do not have data records,then this event is bad
   if NullCount==StaNum  
      eventlist(handles.index).quality=2; 
      Signal=2;
      disp(['All files data in the Fold is Null or NaN:NO.',num2str(handles.index),' Fold:',eventfilename]);
      fid_log=fopen(errorlog,'a+');
      fprintf(fid_log,'\nAll files data in the Fold is Null or NaN:NO.%f\ %s',index,eventfilename);
      fclose(fid_log);           
   end
% sort in GCDkm from small to big
[GCDkm,index]=sort([PAdata.GCDkm]);
PAdata=PAdata(index);

% Plot the waveform file.
function  [ProcessIndex]=PAdataWaveBuide(hObject,eventdata,handles)
global CrossCorrWaveInfo FilterInfo
global PAdata
% global filter cross 
% Sturct creat
PlotWave = struct(CrossCorrWaveInfo);
filter = struct(FilterInfo);

StaNum=size(PAdata,2);
%
PlotWave.BeforeP = str2num(get(handles.BeforePT,'String'));
PlotWave.AfterP = str2num(get(handles.AfterPT,'String'));


RefTravTMin = zeros(1,StaNum);
RefTravTMax = zeros(1,StaNum);
%%%%%%%%%%%%%%%   count the time difference from the the time recorded
% GCDkm/v=sufwave_time+(rcd_time-src_time)
% rcd_Z.DiffT=srctime-rcdtime
% sufwave_time=GCDkm/v+(src_time-rcd_time)

for i = 1:StaNum
     RefTravTMin(i) =  PAdata(i).PT - PlotWave.BeforeP;
     RefTravTMax(i) = PlotWave.AfterP + PAdata(i).PT;
end

if (min(RefTravTMin) <= 0)
    save('RefTimeMin');   
    SaveEvtResult_button_Callback(hObject,eventdata, handles);
    error(['RefTravTMin is negative! NO.',num2str(handles.index),' Fold:',PAdata(1).eventfile]);   
end
StaNum=size(PAdata,2);
time_diff=zeros(1,StaNum);
for loni=1:StaNum
    time_diff(loni)=PAdata(loni).TimeL - RefTravTMax(loni);
    if time_diff(loni)<0
       save('RefTimeMax');
       SaveEvtResult_button_Callback(hObject, eventdata,handles);
       error('The P time window is wrong,the large side of window larger than the record length!');    
    end        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ProcessIndex=1;
% main part 
if ProcessIndex == 1
% preparetion    
   % set filter sample frequency to be the data sampling frequency
   filter.SampleF = PAdata(1).SampleF;
   filter.SampleT = 1/filter.SampleF;
   % make the time window taper
   TaperTime =PlotWave.BeforeP/2 ;  
   TaperNum = round(TaperTime*PAdata(1).SampleF);
   [CrossDeltaT, WaveIndex] = max(RefTravTMax - RefTravTMin);
   % here change cross.StartNum EndNum SENum to PAdata for I dont want
   % cross has to many for 
   for loni=1:StaNum
       PAdata(loni).StartNum = ceil(RefTravTMin(loni)/filter.SampleT);
       PAdata(loni).EndNum = floor(RefTravTMax(loni)/filter.SampleT);
       PAdata(loni).SENum = PAdata(loni).EndNum - PAdata(loni).StartNum + 1;
   end
   PlotWave.PointNum = max([PAdata.SENum]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract the effective surface waves
   % cut-off period sec
   cutoffT=40;
   MaxFilterLength = round(max(512*filter.SampleF, 5*cutoffT*filter.SampleF));
   MaxHalfFilterNum =  floor(MaxFilterLength/2);
   % sort for PAdata and wave in GCDkm
   [GCDkm,index]=sort([PAdata.GCDkm]);
    PAdata=PAdata(index);
    
   for loni=1:StaNum
       PAdata(loni).wave = zeros(1, PlotWave.PointNum + MaxHalfFilterNum); % for the filter
   end
   % set taper window cell 
   window=cell(1,StaNum);
   for loni=1:StaNum
       window{loni}=ones(PAdata(loni).SENum,1);
       window{loni}(1:TaperNum)=sin(0.5*pi*(1:TaperNum)/TaperNum);
       window{loni}((PAdata(loni).SENum-TaperNum+1):PAdata(loni).SENum)=window{loni}(TaperNum:-1:1);
   end
   for loni=1:StaNum
       PAdata(loni).wave(1:PAdata(loni).SENum) = PAdata(loni).OrigWave(PAdata(loni).StartNum:PAdata(loni).EndNum).*window{loni};
   end
end

function LoadFileList(handles,foldname,errorlog)
% 
global EventlistDirectory PAdata eventlist
StaNum=size(PAdata,2);
FileFold=fullfile(EventlistDirectory,foldname);
FileName=[foldname,'_BAD_FileList.txt'];
SavedFile=fullfile(FileFold,FileName);
if exist(SavedFile,'file') 
    if eventlist(handles.index).quality==0
        word=['Handled but  do not on the list! The NO.',num2str(handles.index),' Fold:',eventlist(handles.index).foldname];
        warning(word);
    end
end
if exist(SavedFile,'file')
   count=0;
   SaveFid=fopen(SavedFile,'r');
   while ~feof(SaveFid)
         temp = fscanf(SaveFid, '%s', 1);         
         if ~strcmp(temp,'')
             count=count+1;
             SacList{count}=temp;
         end
   end  
   fclose(SaveFid);
   % assignment
   if count>0
      for loni = 1:count
          for lonj=1:StaNum
              if strcmpi(SacList{loni},PAdata(lonj).fname)
                 PAdata(lonj).quality=2;
              end
          end
      end
   end    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotwaveformAxe1(handles)
global PAdata 
axes1=handles.axes1;
PointNum = max([PAdata.SENum]);
StaNum=size(PAdata,2);
% slider settting
Tlength=handles.Tlength; 
set(handles.SliderPage,'string','1');
if StaNum<=Tlength
    set(handles.PlotSlider,'Max',2);
    set(handles.PlotSlider,'SliderStep',[0. 0.]);
else
    MaxSlider=ceil(StaNum/Tlength);
    set(handles.PlotSlider,'Max',MaxSlider);
    set(handles.PlotSlider,'SliderStep',[1. 1.]/(MaxSlider-1));   
end
Plot_factor=1;   
TraceGap=Plot_factor/2;
hold(axes1,'off');

for loni=1:StaNum
    PlotX=(0:PointNum-1)*PAdata(loni).SampleT;   
    if PAdata(loni).quality==1
      Plot_handle=plot(axes1, PlotX ,TraceGap + PAdata(loni).wave(1:PointNum)/max(PAdata(loni).wave(1:PointNum))*(Plot_factor/2),'k','linewidth',1.4);        
      text('String',PAdata(loni).staname,'Position',[PlotX(end)+1,TraceGap],'horizontalAlignment', 'left','fontsize',9,'fontweight','bold','fontname','Times New Roman','Interpreter','tex','color','k');   
    elseif PAdata(loni).quality==2
      Plot_handle=plot(axes1, PlotX ,TraceGap + PAdata(loni).wave(1:PointNum)/max(PAdata(loni).wave(1:PointNum))*(Plot_factor/2),'c','linewidth',2);   
       text('String',PAdata(loni).staname,'Position',[PlotX(end)+1,TraceGap],'horizontalAlignment', 'left','fontsize',9,'fontweight','bold','fontname','Times New Roman','Interpreter','tex','color','r');   
    end
      PAdata(loni).handle=Plot_handle;      
   
    hold(axes1,'on');
    clear Plot_handle
    TraceGap=TraceGap+Plot_factor;         
end
AxesPlotSetting(axes1);
ylim(axes1,[0 Tlength]);
xlim(axes1,[0 (PointNum-1)*max([PAdata.SampleT])]);
   
% Plot waveform mode
function Plotwaveform_mode(hObject,handles)
global PAdata
fontsize=10;fontweight='bold';fontname='Times New Roman';
PointNum = max([PAdata.SENum]);
StaNum=size(PAdata,2);
% epicenter distance mode plot
Tlength=PAdata(end).GCDdeg-PAdata(1).GCDdeg;
Plot_factor=Tlength/StaNum;   
handles.fig3=figure('Position',[40 46 1600 900],'PaperPosition',[0.1 0.1 10.8 6.075], 'PaperSize',[10 6.275]);
guidata(hObject,handles); 
hold on;
if isempty(PAdata(1).Fwave)
   for loni=1:StaNum
       PlotX=(0:PointNum-1)*PAdata(loni).SampleT;
       if PAdata(loni).quality==1
          plot( PlotX ,PAdata(loni).GCDdeg+ PAdata(loni).wave(1:PointNum)/max(PAdata(loni).wave(1:PointNum))*(Plot_factor/2),'k','linewidth',1.3);   
       elseif PAdata(loni).quality==2
          plot( PlotX ,PAdata(loni).GCDdeg+ PAdata(loni).wave(1:PointNum)/max(PAdata(loni).wave(1:PointNum))*(Plot_factor/2),'c','linewidth',1.8);
       end
    text('String',PAdata(loni).staname,'Position',[PlotX(end)+0.4,PAdata(loni).GCDdeg],'horizontalAlignment', 'left','fontsize',8,'fontweight',fontweight,'fontname',fontname,'Interpreter','tex');      
   end
else
   for loni=1:StaNum
       PlotX=(0:PointNum-1)*PAdata(loni).SampleT;
       if PAdata(loni).quality==1
          plot( PlotX ,PAdata(loni).GCDdeg+ PAdata(loni).Fwave(1:PointNum)/max(PAdata(loni).wave(1:PointNum))*(Plot_factor/2),'k','linewidth',1.3);   
       elseif PAdata(loni).quality==2
          plot( PlotX ,PAdata(loni).GCDdeg+ PAdata(loni).Fwave(1:PointNum)/max(PAdata(loni).wave(1:PointNum))*(Plot_factor/2),'c','linewidth',1.5);
       end
    text('String',PAdata(loni).staname,'Position',[PlotX(end)+0.3,PAdata(loni).GCDdeg],'horizontalAlignment', 'left','fontsize',8,'fontweight',fontweight,'fontname',fontname,'Interpreter','tex');      
   end    
end
current_axes=get(handles.fig3,'CurrentAxes');
AxesPlotSetting(current_axes);
ylabel(current_axes,'Ecpcenter Distance (Deg)');
title(current_axes,[PAdata(1).eventfile(1:8),'  Waveforms in Epicenter Distance'],'fontsize',fontsize,'fontweight',fontweight,'fontname',fontname);
ylim(current_axes,[PAdata(1).GCDdeg-Plot_factor PAdata(end).GCDdeg+Plot_factor]); 
xlim(current_axes,[0 (PointNum-1)*max([PAdata.SampleT])]);
    
% PlotSetting
function AxesPlotSetting(current_axes)

fontsize=10;fontweight='bold';fontname='Times New Roman';
ylabel(current_axes,'Station Trace ','fontsize',fontsize,'fontweight',fontweight,'fontname',fontname);
xlabel(current_axes,'Time (s)','FontSize',fontsize,'FontWeight',fontweight,'fontname',fontname);
xlim(current_axes,[0 3000]); 
% % gca performance
% daspect([50,1,1])
set(current_axes,'linewidth',1.1);
set(current_axes,'fontsize',fontsize,'fontweight',fontweight,'fontname',fontname);

% Updata the message information
function UpdateInform(sta,source,eventlist_filename,handles)
global PAdata eventlist
StaNum=size(PAdata,2);
Ave_epicenter=sum([sta.GCDdeg])/size(sta,2);
Ave_epicenter=(round(Ave_epicenter*1e2))./1e2;
Mag=source.Mag;
set(handles.StaNum,'String',num2str(StaNum));
set(handles.MagEdit,'String',num2str(Mag));
set(handles.DegEdit,'String',num2str(Ave_epicenter));
set(handles.DataEdit,'String',eventlist_filename(3:8));
set(handles.SrcDepthEdit,'String',num2str(sta(1).SrcDepth));
set(handles.SrcLonEdit,'String',num2str(source.Lon))
set(handles.SrcLatEdit,'String',num2str(source.Lat));

index=handles.index;
if eventlist(index).quality==0
    set(handles.Status,'String','New');
elseif eventlist(index).quality==1
    set(handles.Status,'String','Good');
elseif eventlist(index).quality==2
    set(handles.Status,'String','Bad');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Cursor Mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CursorMode(figure)
global dcm_obj
dcm_obj = datacursormode(figure);
set(dcm_obj,'DisplayStyle','datatip',...
    'SnapToDataVertex','off','Enable','on');

function ExitCursorMode()
global dcm_obj
dcm_obj = datacursormode(gcf);
set(dcm_obj,'DisplayStyle','datatip',...
    'SnapToDataVertex','off','Enable','off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Fliter Interface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                         
% --- Executes on button press in Fliter_button.
function Fliter_button_Callback(hObject, eventdata, handles)
% hObject    handle to Fliter_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Tip: the PAdata.wave is cutted by PlotWave function 
FilterProcess(handles);
PlotFilterAxe1(handles);

function FilterProcess(handles)
global PAdata
period1=str2num(get(handles.fq1_edit,'String'));
period2=str2num(get(handles.fq2_edit,'String'));
StaNum=size(PAdata,2);
% dt=1;
for loni=1:StaNum
    PAdata(loni).Fwave=butter_filtfilt(PAdata(loni).wave, PAdata(loni).SampleT,period1,period2);
end

function filter_wave=butter_filtfilt(wave,dt,period1,period2)
butter_num=3;%
if (2 * 1 / period2 * dt)>1 ||(2 * 1 / period1 * dt)>1
   error('Error in Filter! Check your smaple freq and picked filter periods!');
else
  [b,a]=butter(butter_num,[(2 * 1 / period2 * dt),(2 * 1 / period1 * dt)]);
  filter_wave=filtfilt(b,a,wave);
end

function PlotFilterAxe1(handles)
global PAdata 
axes1=handles.axes1;
PointNum = max([PAdata.SENum]);
StaNum=size(PAdata,2);
% slider settting
Tlength=handles.Tlength; 
set(handles.SliderPage,'string','1');
if StaNum<=Tlength
    set(handles.PlotSlider,'Max',2);
    set(handles.PlotSlider,'SliderStep',[0. 0.]);
else
    MaxSlider=ceil(StaNum/Tlength);
    set(handles.PlotSlider,'Max',MaxSlider);
    set(handles.PlotSlider,'SliderStep',[1. 1.]/(MaxSlider-1));   
end
Plot_factor=1;   
TraceGap=Plot_factor/2;
hold(axes1,'off');
for loni=1:StaNum
    PlotX=(0:PointNum-1)*PAdata(loni).SampleT; 
    if PAdata(loni).quality==1
       Plot_handle=plot(axes1, PlotX ,TraceGap + PAdata(loni).Fwave(1:PointNum)/max(PAdata(loni).Fwave(1:PointNum))*(Plot_factor/2),'k','linewidth',1.4);     
       text('String',PAdata(loni).staname,'Position',[PlotX(end)+1,TraceGap],'horizontalAlignment', 'left','fontsize',9,'fontweight','bold','fontname','Times New Roman','Interpreter','tex','color','k');   
    elseif PAdata(loni).quality==2
       Plot_handle=plot(axes1, PlotX ,TraceGap + PAdata(loni).Fwave(1:PointNum)/max(PAdata(loni).Fwave(1:PointNum))*(Plot_factor/2),'c','linewidth',2);
       text('String',PAdata(loni).staname,'Position',[PlotX(end)+1,TraceGap],'horizontalAlignment', 'left','fontsize',9,'fontweight','bold','fontname','Times New Roman','Interpreter','tex','color','r');   
    end
    PAdata(loni).handle=Plot_handle;      

    hold(axes1,'on');
    clear Plot_handle
    TraceGap=TraceGap+Plot_factor;         
end
AxesPlotSetting(axes1);
ylim(axes1,[0 Tlength]);
xlim(axes1,[0 (PointNum-1)*max([PAdata.SampleT])]);


% --- Executes on button press in FilterAuto.
function FilterAuto_Callback(hObject, eventdata, handles)
% hObject    handle to FilterAuto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FilterAuto
if get(handles.FilterAuto,'Value')==1
   set(handles.FilterManual,'Value',0);
   word='Auto Plot the filtered results.';
   disp(word);
   set(handles.MsgEdit,'String',word);
end

% --- Executes on button press in FilterManual.
function FilterManual_Callback(hObject, eventdata, handles)
% hObject    handle to FilterManual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of FilterManual
if get(handles.FilterManual,'Value')==1
   set(handles.FilterAuto,'Value',0);
   word='Manual Plot the filtered results.';
   disp(word);
   set(handles.MsgEdit,'String',word);   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Master  Control Inerface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in Prevbutton.
function Prevbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Prevbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Increase
Increase=-1;
uiresume();
ExitCursorMode
% --- Executes on button press in Goonbutton.
function Goonbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Goonbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Increase
Increase=1;
uiresume();
ExitCursorMode

% --- Executes on button press in Exitbutton.
function Exitbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Exitbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Signal=1, continue Signal=2,skip Signal=3 exit
global Signal
uiresume();
Signal=3;

% --- Executes on button press in SaveEvtResult_button.
function SaveEvtResult_button_Callback(hObject, eventdata,handles)
% hObject    handle to SaveEvtResult_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global EvtsFold
global eventlist
gd_filename='Eventlist_good.txt';
bd_filename='Eventlist_bad.txt';
eventNum=size(eventlist,2);

gd_Series=find([eventlist.quality]==1);
bd_Series=find([eventlist.quality]==2);
NaN_Series=find([eventlist.quality]==0);
if length(NaN_Series)==eventNum
   word=['Non of the event has been markered! Please do PEAnalysis again!'];
   disp(word);
   set(handles.MsgEdit,'String',word);
end
if length(gd_Series)>0
   gd_fid=fopen(fullfile(EvtsFold,gd_filename),'w');
   for loni=1:length(gd_Series)
       fprintf(gd_fid,'%s\n',eventlist(gd_Series(loni)).foldname);
   end
   fclose(gd_fid);
   word=['Good eventlist and bad eventlist files have been rewritten! \n OutPath:',EvtsFold];
   disp(word);
   set(handles.MsgEdit,'String',word);
end
% removed event 
if length(bd_Series)>0
   bd_fid=fopen(fullfile(EvtsFold,bd_filename),'w');
   for loni=1:length(bd_Series)
       fprintf(bd_fid,'%s\n',eventlist(bd_Series(loni)).foldname);
   end
   fclose(bd_fid);
end

% --- Executes on button press in Extract.
function Extract_Callback(hObject, eventdata, handles)
% hObject    handle to Extract (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Extract the good event from goodeventlist and the good sac file from
% orignal fold to the output path.
global EvtsFold EventlistDirectory
gd_filename='Eventlist_good.txt';
gdFile=fullfile(EvtsFold,gd_filename);
% Test weather exist the good eventlist or not
if ~exist(gdFile,'file')
   word='The Output Path May Not Set! Or There is No Good Eventlist!';
   warning(word);
else
   OutPath=fullfile(EvtsFold,'OutFold');
   if ~exist(OutPath,'dir')
       mkdir(OutPath);
   end
   % first level :fold
    [GdFold]=textread(gdFile,'%s');
    if isempty(GdFold)
        warning('The Good Eventlist is Null!');
    else
        disp(['OutPut Path:',EvtsFold]);
        for loni=1:length(GdFold)
              EventFold=fullfile(EventlistDirectory,GdFold{loni});
              FileListPath=fullfile(EventFold,[GdFold{loni},'_FileList.txt']);
              word=['The Event Fold in Good List does not have Filelist Or List is Null:',GdFold{loni}];
              if exist(FileListPath,'file')
                  FileList=textread(FileListPath,'%s');
                  if isempty(FileList)
                      warning(word);
                  else
                      % copy file to output path
                      OutEventFold=fullfile(OutPath,GdFold{loni});
                      if ~exist(OutEventFold,'dir')
                          mkdir(OutEventFold);
                      end
                      for lonj=1:length(FileList)
                            SacName=FileList{lonj};
                            SacFile=fullfile(EventFold,SacName);
                            copyfile(SacFile,OutEventFold);
                      end
                      copyfile(FileListPath,OutEventFold);
                      disp(['Copy Good EventFlod:NO.',num2str(loni),'Fold:',GdFold{loni}]);
                  end
              else                 
                  warning(word);                
              end
        end
    end
end
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Secondary  Control Inerface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in PlotED_button.
function PlotED_button_Callback(hObject, eventdata, handles)
% hObject    handle to PlotED_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 Plotwaveform_mode(hObject,handles);

% --- Executes on button press in Marker_button.
function Marker_button_Callback(hObject, eventdata, handles)
% hObject    handle to Marker_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global PAdata
global dcm_obj
TS_linewidth=2;
c_info = getCursorInfo(dcm_obj);
set(c_info.Target,'LineWidth',TS_linewidth);
set(c_info.Target,'color','c');
curveNum=length(PAdata);
for loni=1:curveNum
    if c_info.Target==PAdata(loni).handle;
       PAdata(loni).quality=2;
    end
end

% --- Executes on button press in Unmark_button.
function Unmark_button_Callback(hObject, eventdata, handles)
% hObject    handle to Unmark_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global PAdata
global dcm_obj
TS_linewidth=1.5;
c_info = getCursorInfo(dcm_obj);
set(c_info.Target,'LineWidth',TS_linewidth);
set(c_info.Target,'color','b');
curveNum=length(PAdata);
for loni=1:curveNum
    if c_info.Target==PAdata(loni).handle;
       PAdata(loni).quality=1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Event Assessment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in Remove_buttion.
function Remove_buttion_Callback(hObject, eventdata, handles)
% hObject    handle to Remove_buttion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global eventlist
index=handles.index;
eventlist(index).quality=2;
word=['The NO.',num2str(handles.index),' Fold:',eventlist(index).foldname,' are removed!'];
set(handles.MsgEdit,'String',word);
disp(word);
uiresume();

% --- Executes on button press in mk_newInfFile_button.
function mk_newInfFile_button_Callback(hObject, eventdata, handles)
% hObject    handle to mk_newInfFile_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Tip: Here the Bad file saved only the handled ones and some bad sac file
% may not read in. So if you want batch the files, please use the Good file
% list.
global eventlist PAdata
global EventlistDirectory

index=handles.index;
eventlist(index).quality =1;
% path
filefold=fullfile(EventlistDirectory,eventlist(index).foldname);
GoodFileName=[eventlist(index).foldname,'_FileList.txt'];
GoodFilePath=fullfile(filefold,GoodFileName);
BadFileName=[eventlist(index).foldname,'_BAD_FileList.txt'];
BadFilePath=fullfile(filefold,BadFileName);
% write 2004 File.txt
GoodFileFig=fopen(GoodFilePath,'w');
BadFileFig=fopen(BadFilePath,'w');
for loni=1:size(PAdata,2)
    if PAdata(loni).quality==1
       fprintf(GoodFileFig,'%s\n',PAdata(loni).fname);
    elseif PAdata(loni).quality==2
       fprintf(BadFileFig,'%s\n',PAdata(loni).fname);
    end
end
fclose(GoodFileFig);
fclose(BadFileFig);




word=['New event has saved NO.',num2str(index),'Fold:',eventlist(index).foldname];
disp(word);
set(handles.MsgEdit,'String',word);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Input & output path Callback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default path for convenience
function DefaultPath
global EvtsFold 
v = version;
disp(['Your current Matlab version ' version]);
if str2double(v(1:4)) < 7.12
   [MPath,junk,junk,junk]=fileparts(mfilename('fullpath'));
else
   [MPath,junk,junk]=fileparts(mfilename('fullpath'));
end
% set the default output paths
OutFold=fullfile(MPath,'OutPut');
if ~exist(OutFold,'dir')
   mkdir(OutFold);
end
EvtsFold =OutFold ;
% add the support path
disp(['Current Path:',MPath]);
addpath(genpath(fullfile(MPath,'SupportLib')));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in eventlist_button.
function eventlist_button_Callback(hObject, eventdata, handles)
% hObject    handle to eventlist_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global EventlistDirectory eventlist 
% eventlist: foldname quality inf eventfile
% global g_SAC_ASCIndex

if EventlistDirectory == 0
   EventlistDirectory = pwd;
end
[pfile1, pname1, index] = uigetfile({'*.txt';'*.dat';'*.*'},'Open file containing all event fold names', EventlistDirectory);
files =[pname1,pfile1];
EventlistDirectory = pname1;
%%%%%%%%%%%%%%%%%%%%%%%%% 
if pname1==0
   word='File containing all event fold names read Unsuccseefully !';
   set(handles.MsgEdit,'String',word);
   warning(word);
else
%%% read the file containing event fold
   fname = fopen(files);
   loni=1;
   while ~feof(fname)
         eventlist(loni).foldname=fgetl(fname);
         eventlist(loni).quality=0; % 0:undetermind; 2:to be removed 1:good, to be saved
         loni=loni+1;
   end
   fclose(fname);
   set(handles.MsgEdit,'String','File containing all event fold names read succseefully !');
end 

% --- Executes on button press in EvntsPath_button.
function EvntsPath_button_Callback(hObject, eventdata, handles)
% hObject    handle to EvntsPath_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global EvtsFold
global eventlist
evtNum=size(eventlist,2);

EvtsFold = uigetdir(EvtsFold, 'Select folder to save and load the output files:');
if EvtsFold==0
   EvtsFold= fullfile(pwd,'OutPut');
end
set(handles.MsgEdit,'String',['Output folder is: ',EvtsFold]);
% search the good.txt bad.txt
% gd_filename='Eventlist_good.txt';
% bd_filename='Eventlist_bad.txt';
gd_keyword='*good.txt';
gd_dir=dir(fullfile(EvtsFold,gd_keyword)); 
% for good eventlist
if ~isempty(gd_dir) && (~isempty(eventlist))
    gd_fname=gd_dir.name;
    gd_file=fullfile(EvtsFold,gd_fname);
    disp('The existing Good event list file has benn read!');
   % load list
   count=0;
   gd_fid=fopen(gd_file,'r');
   while ~feof(gd_fid)
         temp = fscanf(gd_fid, '%s', 1);
         count=count+1;
         if ~strcmp(temp,'')
             gd_list{count}=temp;
         end
   end  
   fclose(gd_fid);
   % assignment
   gd_evtNum=size(gd_list,2);
   if gd_evtNum>0
      for loni = 1:gd_evtNum
          for lonj=1:evtNum
              if strcmpi(gd_list{loni},eventlist(lonj).foldname)
                 eventlist(lonj).quality=1;
              end
          end
      end
   end    
end

% for bad eventlist
bd_keyword='*bad.txt';
bd_dir=dir(fullfile(EvtsFold,bd_keyword));
if ~isempty(bd_dir) && (~isempty(eventlist))
   bd_fname=bd_dir.name;
   bd_file=fullfile(EvtsFold,bd_fname);
   count=0;
   disp('The existing Bad event list file has benn read!');
   bd_fid=fopen(bd_file,'r');
   while ~feof(bd_fid)
         temp = fscanf(bd_fid, '%s', 1);
         count=count+1;
         if ~strcmp(temp,'')
             bd_list{count}=temp;
         end
   end  
   fclose(bd_fid);
   bd_evtNum=size(bd_list,2);
   if bd_evtNum>0
      for loni = 1:bd_evtNum
          for lonj=1:evtNum
              if strcmpi(bd_list{loni},eventlist(lonj).foldname)
                 eventlist(lonj).quality=2;
              end
          end
      end
   end    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       MsEdit Callback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MsgEdit_Callback(hObject, eventdata, handles)
% hObject    handle to MsgEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of MsgEdit as text
%        str2double(get(hObject,'String')) returns contents of MsgEdit as a double
% --- Executes during object creation, after setting all properties.

function MsgEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MsgEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Unicontrol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on slider movement.
function PlotSlider_Callback(hObject, eventdata, handles)
% hObject    handle to PlotSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
SliderIndex=get(handles.PlotSlider,'Value');
SliderIndex=round(SliderIndex);
set(handles.SliderPage,'string',num2str(SliderIndex));
Tlength=handles.Tlength;
ylimRange=[0 Tlength]+Tlength*(SliderIndex-1);
axes1=handles.axes1;
ylim(axes1,ylimRange);
% --- Executes during object creation, after setting all properties.
function PlotSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PlotSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in NoFliter.
function NoFliter_Callback(hObject, eventdata, handles)
% hObject    handle to NoFliter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PlotwaveformAxe1(handles);

% --- Executes on button press in PlotStation.
function PlotStation_Callback(hObject, eventdata, handles)
% hObject    handle to PlotStation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% push the botton after load the wave file!
global PAdata
load topomap 
fontsize=10;fontweight='bold';fontname='Times New Roman';
fig4=figure('Position',[40 46 1600 900],'PaperPosition',[0.1 0.1 10.8 6.075], 'PaperSize',[7 6]);
axes=gca;
x=linspace(0,360,1080);
y=linspace(-90,90,2160);
[cmap clim] = demcmap(topomap);
imagesc(x,y,topomap,clim);
ColBarH=colorbar('vert');
set(get(ColBarH,'title'),'string','Topo (km)','fontsize',fontsize,'fontweight',fontweight,'fontname',fontname);
hold(axes,'on');
colormap(cmap);axis image; grid on; axis on;
load coast 
kk = find(long < 0);
long(kk) = 360 - abs(long(kk));
ii = find( abs(long) < 0.5);
long(ii) = NaN;
lat(ii) = NaN;
plot(axes,long,lat,'k', 'LineWidth',2);
set(axes,'ydir','normal');
set(axes,'FontSize',16,'FontWeight','bold');
axis equal

stanum = size(PAdata,2);
Lon=zeros(1,stanum);
Lat=zeros(1,stanum);
for i = 1:stanum
    Lon(i) = PAdata(i).Lon;
    Lat(i) = PAdata(i).Lat;
    plot(axes,  Lon(i), Lat(i), 'k^', 'MarkerSize',6, 'MarkerFaceColor','k');
end
ExtraLon = 0.1*(max(Lon)-min(Lon));
ExtraLat = 0.1*(max(Lat)-min(Lat));
ExtraLon = min(ExtraLon, 1);
ExtraLat = min(ExtraLat, 1);
xlim(axes,[min(Lon)-ExtraLon max(Lon)+ExtraLon]);
ylim(axes,[min(Lat)-ExtraLat max(Lat)+ExtraLat]);
set(axes, 'XTickMode','auto','YTickMode','auto');
title(axes,'Station Map','fontsize',fontsize,'fontweight',fontweight,'fontname',13);
ylabel(axes,'Latitude ','fontsize',fontsize,'fontweight',fontweight,'fontname',fontname);
xlabel(axes,'Longitude','FontSize',fontsize,'FontWeight',fontweight,'fontname',fontname);
set(axes,'linewidth',1.1);
set(axes,'fontsize',fontsize,'fontweight',fontweight,'fontname',fontname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Creat Fcn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
