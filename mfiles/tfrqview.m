function tfrqview(tfr,sig,t,method,p1,p2,p3,p4,p5);
%TFRQVIEW Quick visualization of time-frequency representations.
%       TFRQVIEW(TFR,SIG,T,METHOD,P1,P2,P3,P4,P5) allows a quick 
%       visualization of a time-frequency representation.
%
%       TFR     : time-frequency representation (MxN).
%       SIG     : signal in time. If unavailable, put sig=[] as input
%                 parameter.                    (default : []).
%       T       : time instants                 (default : 1:N).
%       METHOD  : name of chosen representation (default : 'TYPE1').
%                 See the TFR* files for authorized names.
%                 TYPE1 : the representation TFR goes in normalized
%                       frequency from -0.5 to 0.5 ; 
%                 TYPE2 : the representation TFR goes in normalized
%                       frequency from 0 to 0.5. 
%       P1...P5 : optional parameters of the representation : run the 
%                 file TFRPARAM(METHOD) to know the meaning of P1..P5 
%                 for your method. 
%
%	When you use the 'save' option in the main menu, you save all your
%	variables as well as two strings, TfrQView and TfrView, in a mat 
%	file. If you load this file and do eval(TfrQView), you will restart
%	the display session under tfrqview ; if you do eval(TfrView), you
%	will obtain the exact layout of the screen you had when clicking on 
%	the 'save' button. 
%   
%       Example : 
%        sig=fmsin(128); tfr=tfrwv(sig);
%        tfrqview(tfr,sig,1:128,'tfrwv');
%
%       See also TFRVIEW, TFRSAVE, TFRPARAM.

%       F. Auger, September 1994, July 1995 
%       O. Lemoine, Oct 1995, May-July 1996.
%       F. Auger, May 1998.
%       Copyright (c) 1996 by CNRS (France).
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

comp=computer;  %  so as to know the running computer
MatlabVersion=version; MatlabVersion=str2num(MatlabVersion(1));

% Tests on the input arguments
if nargin<1,
 error('At least one parameter required'); % at least the tfr
end
[tfrrow,tfrcol]=size(tfr);
if nargin==1,
 sig=[]; t=1:tfrcol; method='type1'; % empty signal
elseif nargin==2,
 t=1:tfrcol; method='type1';
elseif nargin==3,
 method='type1';
end
[trow,tcol] = size(t);
if (trow~=1),
 error('T must only have one row'); % t must be a row vector
end;
if (tfrcol~=tcol),
 error('T must have as much elements as tfr has columns');
end;

[Nsig,Ncol]=size(sig);
if Ncol>2,
 error('SIG must have one or two columns');
end

% Computation of Nf2, the number of interresting points in frequency
method=upper(method);
if istfr2(method),
 Nf2=tfrrow;
elseif istfr1(method),
 Nf2=tfrrow/2; 
else
 error('Unknown representation. Use type1 or type2');
end;

% Computation of freq (vector of frequency samples) 
if istfraff(method),
 freq=eval(['p',num2str(nargin-4)]);   	% last input argument is freqs.
else
 freq=(0.5*(0:Nf2-1)/Nf2);
end

% Initialization of the variables
if exist('options.mat'),
 load options
 colormap(SavedColorMap);
else
 threshold=5.0;  % visualization threshold
 linlogtfr=0;    % tfr visualization scale : 0 for linear 1 for logarithmic
 linlogspec=1;   % spectrum visualization scale
 sigenveloppe=0; % signal enveloppe visualization

 levelnumb=64;   % number of levels in the contour plot
 colmap=1;       % colormap index
 
 display=2;      % display index
 
 isgridsig=0;    % grid on signal
 isgridspec=0;   % grid on spectrum
 isgridtfr=0;    % grid on tfr
 
 issig=0;        % display signal
 isspec=0;       % display spectrum
 iscolorbar=0;   % display colorbar

 fs=1.0;         % sampling frequency (Hz)
 fmin=0.0;       % smallest displayed frequency 
 fmax=0.5*fs;    % highest displayed frequency
end;


% Test of analycity
if ~isempty(sig),
 for k=1:Ncol,
  % spec(:,k)=abs(fft(sig(min(t):max(t),k))).^2; Nsp=length(spec);
  % modifications :  F. Auger (fog), 30/11/97
  Lt_fog=max(t)-min(t)+1;   
  Nb_tranches_fog = floor(Lt_fog/tfrrow);
  % fprintf('%f \n',Nb_tranches_fog);
  spec(:,k)=zeros(tfrrow,1);
  for Num_tranche_fog=0:Nb_tranches_fog-1,
   % fprintf('%f \n',Num_tranche_fog);
   spec(:,k)=spec(:,k)+abs(fft(sig(min(t)+tfrrow*Num_tranche_fog+(0:tfrrow-1),k))).^2;
  end;

  if (Lt_fog>Nb_tranches_fog*tfrrow),
   spectre_fog=fft(sig(min(t)+tfrrow*Nb_tranches_fog:max(t),k),tfrrow);
   spectre_fog=spectre_fog(:);
   spec(:,k)=spec(:,k)+abs(spectre_fog).^2;
  end;

  
  % spec1=sum(spec(1:tfrrow/2,k));
  % spec2=sum(spec(tfrrow/2+1:Nsp,k));
  spec1=sum(spec(1:tfrrow/2,k));
  spec2=sum(spec(tfrrow/2+1:tfrrow,k));
  
  if spec2>spec1/10,
   disp('Be careful : the signal is not analytic!');
  end
 end
end

% Test of reality
if (Ncol==2 & ~isreal(tfr)),
 disp('Cross distribution. As the result is complex, we display the real part.');
 tfr=real(tfr);
end

ChoiceDisplay     =  1; % All the possible values of the choice variable
ChoiceLayout      =  2;
ChoiceSampling    =  3;
ChoiceFreqBounds  =  4;
ChoiceThreshold   =  5;
ChoiceLinlog      =  6;
ChoiceRedraw      =  7;
ChoiceNewFigure   =  8;
ChoiceSaveResults =  9;
ChoiceSaveOptions = 10;
ChoicePrint       = 11;
ChoiceClose       = 12;

CallTfrView = 1; % 1 to call tfrview, 0 not to do it
RefreshFigure=1; % 1 to refresh figure every time, 0 to freeze

choice=ChoiceSampling;
while choice~=ChoiceClose,                       % while not close
 if RefreshFigure & CallTfrView,                 % Call to tfrview
  linlog=linlogtfr+2*linlogspec+4*sigenveloppe;
  isgrid=isgridsig+2*isgridspec+4*isgridtfr;
  layout=issig+isspec*2+iscolorbar*4+1;
  param = [display, linlog, threshold, levelnumb, Nf2, layout,...
           fs, isgrid, fmin, fmax];
  if (nargin<=4),
   tfrview(tfr,sig,t,method,param);
  elseif (nargin==5),
   tfrview(tfr,sig,t,method,param,p1);
  elseif (nargin==6),
   tfrview(tfr,sig,t,method,param,p1,p2);
  elseif (nargin==7),
   tfrview(tfr,sig,t,method,param,p1,p2,p3);
  elseif (nargin==8),
   tfrview(tfr,sig,t,method,param,p1,p2,p3,p4);
  elseif (nargin==9),
   tfrview(tfr,sig,t,method,param,p1,p2,p3,p4,p5);
  end;
 end;

 if (linlogtfr==0),                              % Lin/log scale of the tfr
  linlogstr='Change to a logarithmic scale';
 else
  linlogstr='Change to a linear scale';
 end;

 if (RefreshFigure==1),
  redrawstr='Don''t redraw yet';
 else
  redrawstr='Redraw now';
 end;

 % Main menu
 choice=menu ('TFRQVIEW MENU :',...
              'Change the display mode',...       % ChoiceDisplay
              'Change the display layout',...     % ChoiceLayout
              'Change the sampling frequency',... % ChoiceSampling
              'Change the frequency bounds',...   % ChoiceFreqBounds
              'Change the threshold',...          % ChoiceThreshold
              linlogstr,...                       % ChoiceLinlog
              redrawstr,...                       % ChoiceRedraw
              'New figure',...                    % ChoiceNewFigure
              'Save results',...                  % ChoiceSaveResults
              'Save options',...                  % ChoiceSaveOptions
              'Print',...                         % ChoicePrint
              'Close');                           % ChoiceClose
 
 if (choice==ChoiceDisplay),                      % Change the display mode

  OldDisplay=display;
  display=menu('DISPLAY MODE :',...
               'contour',...                      % 1
               'imagesc',...                      % 2
               'pcolor',...                       % 3
               'surf',...                         % 4
               'mesh',...                         % 5
               'change the color map',...         % 6
               'change the number of colors or levels',...  % 7
               'cancel');                         % 8

  if (display>=1)&(display<=5),
   CallTfrView=1;
  elseif (display==6),
   if MatlabVersion>=5,
    colmap=menu('COLOR MAP :',...
                'hsv','jet','cool','bone','gray','hot','prism',...
                'pink','colorcube','autumn','winter','spring','summer',...
                'brighten','darken','permute','spin','cancel');
    if     colmap== 1,  colormap(hsv(levelnumb));
    elseif colmap== 2,  colormap(jet(levelnumb));
    elseif colmap== 3,  colormap(cool(levelnumb));
    elseif colmap== 4,  colormap(bone(levelnumb));
    elseif colmap== 5,  colormap(gray(levelnumb));
    elseif colmap== 6,  colormap(hot(levelnumb));
    elseif colmap== 7,  colormap(prism(levelnumb));
    elseif colmap== 8,  colormap(pink(levelnumb));
    elseif colmap== 9,  colormap(colorcube(levelnumb));
    elseif colmap==10,  colormap(autumn(levelnumb));
    elseif colmap==11,  colormap(winter(levelnumb));
    elseif colmap==12,  colormap(spring(levelnumb));
    elseif colmap==13,  colormap(summer(levelnumb));
    elseif colmap==14,  brighten(+0.20);
    elseif colmap==15,  brighten(-0.10);
    elseif colmap==16,  MyMap = colormap; colormap(flipud(MyMap));
    elseif colmap==17,  spinmap;
    end
   else
    colmap=menu('COLOR MAP :',...
                'hsv','jet','cool','bone','gray','hot','prism',...
                'brighten','darken','permute','spin','cancel');
    if     colmap== 1,  colormap(hsv(levelnumb));
    elseif colmap== 2,  colormap(jet(levelnumb));
    elseif colmap== 3,  colormap(cool(levelnumb));
    elseif colmap== 4,  colormap(bone(levelnumb));
    elseif colmap== 5,  colormap(gray(levelnumb));
    elseif colmap== 6,  colormap(hot(levelnumb));
    elseif colmap== 7,  colormap(prism(levelnumb));
    elseif colmap== 8,  brighten(+0.25);
    elseif colmap== 9,  brighten(-0.25);
    elseif colmap==10,  MyMap = colormap; colormap(flipud(MyMap));
    elseif colmap==11,  spinmap;
    end
   end
   display=OldDisplay; CallTfrView=0;

  elseif (display==7),
   fprintf(' Old number of levels: %f\n',levelnumb); levelold=levelnumb;
   levelnumb=input(' New number of levels: ');
   if isempty(levelnumb),
    levelnumb=levelold; CallTfrView=0;
   else
    if levelnumb<levelold,
     CallTfrView=1;
     MyMap = colormap; MyMap=MyMap(1:levelnumb,:); colormap(MyMap);
    elseif levelnumb>levelold,
     CallTfrView=1;
     MyMap = ones(levelnumb, 3); MyMap(1:levelold,:)=colormap; 
     fprintf('warning : the colormap size has been increased by identical vectors\n');
     fprintf('You should redefine the colormap\n');
    else
     CallTfrView=0;
    end
   end
   display=OldDisplay; 

  elseif (display==8),
   display=OldDisplay; CallTfrView=0;
  end;

 elseif (choice==ChoiceLayout),                 % Change the display layout
 
  layout=1;
  if issig==0, 
   SignalStr= 'display signal';
  else
   SignalStr='remove signal';
  end;
 
  if isspec==0, 
   SpectrumStr= 'display spectrum';
  else
   SpectrumStr='remove spectrum';
  end;

  if ~issig & ~isspec,
   if isgridtfr,
    GridStr='Remove the grid';
   else
    GridStr='Add a grid';
   end;
  else
   GridStr='Grids';
  end;
 
  if iscolorbar==0,
   ColorbarStr='display colorbar';
  else
   ColorbarStr='remove colorbar';
  end;
  
  layout=menu('DISPLAY LAYOUT',...
               SignalStr,...
               SpectrumStr,...
               GridStr,...
               ColorbarStr,...
              'cancel');
            
  if layout==1,
   issig=~issig;
   if issig==1, 
    if isempty(sig),
     disp('Impossible action : the signal is unavailable'); issig=0; CallTfrView=0;
    else
     sigenveloppe=menu('SIGNAL REPRESENTATION','signal only','signal with enveloppe')-1;
     CallTfrView=1;
    end;
   else
    isgridsig=0;
   end; 
  elseif layout==2,   
   isspec=~isspec;
   if isspec==1,
    if isempty(sig),
     disp('Impossible action : the signal is unavailable'); isspec=0; CallTfrView=0;
    else    
     linlogspec=menu('FREQUENCY REPRESENTATION','linear scale','log scale')-1;
     CallTfrView=1;
    end;
   else
    isgridspec=0;
   end;

  elseif layout==3,

   if ~issig & ~isspec,	                 % No signal and no spectrum
    isgridtfr=1-isgridtfr; 
    CallTfrView=1;  

   elseif issig & ~isspec,               % A signal, no spectrum 
    Grid=1;
    if ~isgridsig,
     gridsigstr='add a grid on the signal';
    else
     gridsigstr='remove the grid on the signal';
    end

    if ~isgridtfr,
     gridtfrstr='add a grid on the TFR';
    else
     gridtfrstr='remove the grid on the TFR';
    end

    Grid=menu('GRID MENU :',gridsigstr,gridtfrstr,'cancel');
    if Grid==1,
     isgridsig=1-isgridsig; CallTfrView=1;
    elseif Grid==2,
     isgridtfr=1-isgridtfr; CallTfrView=1;
    else
     CallTfrView=0;
    end

   elseif ~issig & isspec,               % No signal, a spectrum
    Grid=1;
    if ~isgridspec,
     gridspestr='add a grid on the spectrum';
    else
     gridspestr='remove the grid on the spectrum';
    end
    if ~isgridtfr,
     gridtfrstr='add a grid on the TFR';
    else
     gridtfrstr='remove the grid on the TFR';
    end
    Grid=menu('GRID MENU :',gridspestr,gridtfrstr,'Close');
    if Grid==1,
     isgridspec=1-isgridspec; CallTfrView=1;
    elseif Grid==2,
     isgridtfr=1-isgridtfr; CallTfrView=1;
    else CallTfrView=0;
    end

  else                                  % A signal and a spectrum
   Grid=1;
   if ~isgridsig,
    gridsigstr='add a grid on the signal';
   else
    gridsigstr='remove the grid on the signal';
   end
   if ~isgridspec,
    gridspestr='add a grid on the spectrum';
   else
    gridspestr='remove the grid on the spectrum';
   end
   if ~isgridtfr,
    gridtfrstr='add a grid on the TFR';
   else
    gridtfrstr='remove the grid on the TFR';
   end
   Grid=menu('GRID MENU :',gridsigstr,gridspestr,gridtfrstr,'cancel');
   if Grid==1,
    isgridsig=1-isgridsig; CallTfrView=1;
   elseif Grid==2,
    isgridspec=1-isgridspec; CallTfrView=1;
   elseif Grid==3,
    isgridtfr=1-isgridtfr; CallTfrView=1;
   else CallTfrView=0;
   end
  end


  elseif layout==4,
   iscolorbar=~iscolorbar; CallTfrView=1;
  elseif layout==5,
   CallTfrView=0;
  end;             
  
 elseif (choice==ChoiceSampling),                   % Change the sampling frequency 

  fprintf(' Old sampling frequency: %f\n',fs); 
  fsold=fs; fs=input(' New sampling frequency: ');
  if isempty(fs),
   fs=fsold; CallTfrView=0;
  else
   CallTfrView=1;
  end; 
 elseif (choice==ChoiceFreqBounds),                 % Change the frequency bounds
  CallTfrView=0;

     fprintf(' Old smallest normalized frequency : %f\n',fmin); fminold=fmin; 
  fmin=input(' New smallest normalized frequency : ');
  if isempty(fmin),
   fmin=fminold; 
  elseif fmin>0.5,
   fprintf('normalized frequency desired ! value unmodified\n');
   fmin=fminold; 
  else
   CallTfrView=1;
  end; 

     fprintf(' Old highest normalized frequency  : %f\n',fmax); fmaxold=fmax; 
  fmax=input(' New highest normalized frequency  : ');
  if isempty(fmax),
   fmax=fmaxold; 
  elseif fmax>0.5,
   fprintf('normalized frequency desired ! value unmodified\n');
   fmax=fmaxold;
  else
   CallTfrView=1;
  end; 
  
 elseif (choice==ChoiceThreshold),                  % Change the threshold

  fprintf(' Old threshold: %f\n', threshold); throld=threshold;
  threshold=input(' New threshold: ');
  if isempty(threshold),
   threshold=throld; CallTfrView=0;
  else
   CallTfrView=1;
  end

 elseif (choice==ChoiceLinlog),                     % Change the lin/log scale of tfr

  linlogtfr=1-linlogtfr;

 elseif (choice==ChoiceRedraw),                           % redraw ?
  RefreshFigure=1-RefreshFigure;
  if RefreshFigure==1, CallTfrView=1; end;

 elseif (choice==ChoiceNewFigure),                        % new figure
  figure; CallTfrView=1;

 elseif (choice==ChoiceSaveResults),                      % Save the results

  f=freq*fs;
  Nmethod=length(method);
  if comp(1:2)=='PC',
   DefaultName=[method(4:Nmethod),num2str(Nsig),'.mat'];
   [name,PathWorkDir] = uiputfile(DefaultName, 'Save As');
  else
   DefaultName=[method(4:Nmethod),num2str(Nsig)];
   nameStr=[' Name of the MAT file [',DefaultName,'] : '];
   name=input(nameStr,'s'); 
   while (length(name)>8),
    disp(' The name must have less than 8 characters');
    name=input(nameStr,'s'); 
   end
   if isempty(name),
    name=DefaultName;
   end
   PathWorkDir='';
  end
  linlog=linlogtfr+2*linlogspec+4*sigenveloppe;
  isgrid=isgridsig+2*isgridspec+4*isgridtfr;
  param = [display,linlog,threshold,levelnumb,Nf2,layout,fs,isgrid,fmin,fmax];
  SavedColorMap=colormap;
  if (nargin<=4),
   TfrQView=['colormap(SavedColorMap); tfrqview(tfr,sig,t,method)'];
   TfrView =['clf;colormap(SavedColorMap); tfrview(tfr,sig,t,method,param)'];
   eval(['save ',PathWorkDir,name,...
         ' tfr sig t f fs method param SavedColorMap TfrView TfrQView']);
  elseif (nargin==5),
   TfrQView=['colormap(SavedColorMap); tfrqview(tfr,sig,t,method,p1)'];
   TfrView =['clf; colormap(SavedColorMap); tfrview(tfr,sig,t,method,param,p1)']; 
   eval(['save ',PathWorkDir,name, ...
         ' tfr sig t f fs method param p1 SavedColorMap TfrView TfrQView']);
  elseif (nargin==6),
   TfrQView=['colormap(SavedColorMap); tfrqview(tfr,sig,t,method,p1,p2)'];
   TfrView =['clf; colormap(SavedColorMap); tfrview(tfr,sig,t,method,param,p1,p2)'];
   eval(['save ',PathWorkDir,name,...
         ' tfr sig t f fs method param p1 p2 SavedColorMap TfrView TfrQView']);
  elseif (nargin==7),
   TfrQView=['colormap(SavedColorMap); tfrqview(tfr,sig,t,method,p1,p2,p3)'];
   TfrView =['clf; colormap(SavedColorMap); tfrview(tfr,sig,t,method,param,p1,p2,p3)'];
   eval(['save ',PathWorkDir,name,...
         ' tfr sig t f fs method param p1 p2 p3 SavedColorMap TfrView TfrQView']);
  elseif (nargin==8),
   TfrQView=['colormap(SavedColorMap); tfrqview(tfr,sig,t,method,p1,p2,p3,p4)'];
   TfrView =['clf; colormap(SavedColorMap); tfrview(tfr,sig,t,method,param,p1,p2,p3,p4)'];
   eval(['save ',PathWorkDir,name,...
         ' tfr sig t f fs method param p1 p2 p3 p4 SavedColorMap TfrView TfrQView']);
  elseif (nargin==9),
   TfrQView=['colormap(SavedColorMap); tfrqview(tfr,sig,t,method,p1,p2,p3,p4,p5)'];
   TfrView =['clf; colormap(SavedColorMap); tfrview(tfr,sig,t,method,param,p1,p2,p3,p4,p5)'];
   eval(['save ',PathWorkDir,name,...
         ' tfr sig t f fs method param p1 p2 p3 p4 p5 SavedColorMap TfrView TfrQView']);
  end;
  disp(' ');
  fprintf('The file is saved in the directory %s\n',PathWorkDir);
  fprintf('under the name %s\n',name);

  fprintf('If you want to find again the exact layout of this screen, do\n');
  fprintf('load %s%s; eval(TfrView);\n\n', PathWorkDir,name);
  fprintf('If you want to restart the display session under tfrqview, do\n');
  fprintf('load %s%s; eval(TfrQView);\n',PathWorkDir,name);
  CallTfrView=0;

 elseif (choice==ChoiceSaveOptions),                 % Save options
  SavedColorMap=colormap;
  save options fs fmin fmax threshold linlogtfr linlogspec levelnumb  ...
       display layout colmap SavedColorMap iscolorbar  ...
       isgridsig isgridspec isgridtfr issig sigenveloppe isspec;
  fprintf('\n Options saved\n');
  CallTfrView=0;
 
 elseif (choice==ChoicePrint),	                     % Print the current figure

  Nmethod=length(method);
  TFTBDevice = MENU('Choose a device',...
                    '-deps','-depsc','-deps2','-depsc2','-djpeg','-dtiff','other','cancel');
  if TFTBDevice==1,
   TFTBDeviceName='-deps '  ; TFTBExtension='.eps';
  elseif TFTBDevice==2,
   TFTBDeviceName='-depsc ' ; TFTBExtension='.eps';
  elseif TFTBDevice==3,
   TFTBDeviceName='-deps2 ' ; TFTBExtension='.eps';
  elseif TFTBDevice==4,
   TFTBDeviceName='-depsc2 '; TFTBExtension='.eps';
  elseif TFTBDevice==5,
   TFTBDeviceName='-djpeg ' ; TFTBExtension='.jpg';
  elseif TFTBDevice==6,
   TFTBDeviceName='-dtiff ' ; TFTBExtension='.tif';
  elseif TFTBDevice==7,
   TFTBDeviceName=input('device option : ','s'); ; TFTBDeviceName=[TFTBDeviceName,' '];
   TFTBExtension =input('file extension : ','s'); ;
  end;
  if TFTBDevice~=8,
   if comp(1:2)=='PC',
    DefaultName=[method(4:Nmethod),num2str(Nsig),TFTBExtension];
    [name,PathWorkDir] = uiputfile(DefaultName, 'Save As');
   else
    DefaultName=[method(4:Nmethod),num2str(Nsig),TFTBExtension];
    nameStr=[' file name [',DefaultName,'] : '];
    name=input(nameStr,'s'); 
    while (length(name)>8),
     disp('The name must have less than 8 characters');
     name=input(nameStr,'s'); 
    end
    if isempty(name),
     name=DefaultName;
    end
   end
   % ['print ', TFTBDeviceName, PathWorkDir, name]
   eval(['print ', TFTBDeviceName, PathWorkDir, name]); 
   fprintf(' The file is saved in the directory %s\n',PathWorkDir);
   fprintf('under the name %s\n',name);
  end;
  CallTfrView=0;

 end;

end;

% good bye. I hope that everything happened fine.
