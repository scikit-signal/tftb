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

comp=computer; 

% Tests on the input arguments
if nargin<1,
 error('At least one parameter required');
end
[tfrrow,tfrcol]=size(tfr);
if nargin==1,
 sig=[]; t=1:tfrcol; method='type1';
elseif nargin==2,
 t=1:tfrcol; method='type1';
elseif nargin==3,
 method='type1';
end
[trow,tcol] = size(t);
if (trow~=1),
 error('T must only have one row'); 
end;
if (tfrcol~=tcol),
 error('T must have as much elements as tfr has columns');
end;
method=upper(method);
[Nsig,Ncol]=size(sig);
if Ncol>2,
 error('SIG must have less than two columns');
end

% Computation of Nf2, the number of interresting points in frequency
method=upper(method);
if strcmp(method,'TFRWV'   ) | strcmp(method,'TFRPWV'  ) | ...
   strcmp(method,'TFRSPWV' ) | strcmp(method,'TFRCW'   ) | ...
   strcmp(method,'TFRZAM'  ) | strcmp(method,'TFRBJ'   ) | ...
   strcmp(method,'TFRBUD'  ) | strcmp(method,'TFRGRD'  ) | ...
   strcmp(method,'TFRRSPWV') | strcmp(method,'TFRRPWV' ) | ...
   strcmp(method,'TFRRIDB' ) | strcmp(method,'TFRRIDH' ) | ...
   strcmp(method,'TFRRIDT' ) | strcmp(method,'TFRASPW' ) | ...
   strcmp(method,'TFRDFLA' ) | strcmp(method,'TFRSPAW' ) | ...
   strcmp(method,'TFRRIDBN') | strcmp(method,'TFRUNTER') | ...
   strcmp(method,'TFRBERT' ) | strcmp(method,'TFRSCALO') | ...
   strcmp(method,'TYPE2'   ),
 Nf2=tfrrow;
elseif strcmp(method,'TFRPMH'  )| strcmp(method,'TFRRPMH' )| ...
       strcmp(method,'TFRSP'   )| strcmp(method,'TFRRSP'  )| ...
       strcmp(method,'TFRPPAGE')| strcmp(method,'TFRRPPAG')| ...
       strcmp(method,'TFRMHS'  )| strcmp(method,'TFRRGAB' )| ...
       strcmp(method,'TFRMH'   )| strcmp(method,'TFRMMCE' )| ...
       strcmp(method,'TFRRMSC' )| strcmp(method,'TFRPAGE' )| ...
       strcmp(method,'TFRGABOR')| strcmp(method,'TFRRI'   )| ...
       strcmp(method,'TFRMSC'  )| strcmp(method,'TYPE1'   )| ...
       strcmp(method,'TFRSTFT' ),
 Nf2=tfrrow/2; 
else
 error('Unknown representation. Use type1 or type2');
end;

% Computation of freq (vector of frequency samples) 
if strcmp(method,'TFRASPW' ) | strcmp(method,'TFRSCALO') | ...
   strcmp(method,'TFRDFLA' ) | strcmp(method,'TFRSPAW' ) | ...
   strcmp(method,'TFRUNTER') | strcmp(method,'TFRBERT' ),
 freq=eval(['p',num2str(nargin-4)]);   	% last input argument is freqs.
else
 freq=(0.5*(0:Nf2-1)/Nf2);
end


% Initialization of the current figure
zoom off; clf; 
set(gcf,'Resize','On','NextPlot','Add');

% Initialization of the variables
if exist('options.mat'),
 load options
else
 threshold=5.0;
 linlogtfr=0;
 linlogspec=1;
 levelnumb=6;
 colmap=1;
 display=1;
 layout=1;
 iscolorbar=0; 
 isgridsig=0; 
 isgridspe=0; 
 isgridtfr=0;
 issig=0; 
 isspec=0;
end
choice=4;
Grid=1;
displayold=1;
fs=1.0;

if ((layout>=2 & layout<=4) & isempty(sig)),
 layout=1;
end


% Initialization of the axes 
axcb   = axes('Units','normal','Visible','off','Box','On');
axsig  = axes('Units','normal','Visible','off','Box','On');
axspec = axes('Units','normal','Visible','off','Box','On');
axtfr  = axes('Units','normal','Visible','off','Box','On');
if comp(1:2)=='PC',
 set(axsig ,'fontsize',10); 
 set(axspec,'fontsize',10); 
 set(axtfr ,'fontsize',10); 
end
set(gcf,'UserData',[axcb axsig axspec axtfr]);

% Test of analycity
if ~isempty(sig),
 for k=1:Ncol,
  spec(:,k)=abs(fft(sig(min(t):max(t),k))).^2; Nsp=length(spec);
  spec1=sum(spec(1:round(Nsp/2),k));
  spec2=sum(spec(round(Nsp/2)+1:Nsp,k));
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


while choice~=12,

 axes(axtfr);				% The current axis is the tfr axis

 if (display~=7 & choice~=2 & choice~=3 & choice<8),	% Call to tfrview
  linlog=linlogtfr+2*linlogspec;
  if choice==4,
   access=3;	
  else
   access=1;
  end
  state=issig+2*iscolorbar;
  isgrid=isgridsig+2*isgridspe+4*isgridtfr;
  param = [display,linlog,threshold,levelnumb,Nf2,layout,access,state,fs,isgrid];
  map=colormap;
  if (nargin<=4),
   tfrview(tfr,sig,t,method,param,map);
  elseif (nargin==5),
   tfrview(tfr,sig,t,method,param,map,p1);
  elseif (nargin==6),
   tfrview(tfr,sig,t,method,param,map,p1,p2);
  elseif (nargin==7),
   tfrview(tfr,sig,t,method,param,map,p1,p2,p3);
  elseif (nargin==8),
   tfrview(tfr,sig,t,method,param,map,p1,p2,p3,p4);
  elseif (nargin==9),
   tfrview(tfr,sig,t,method,param,map,p1,p2,p3,p4,p5);
  end;
 elseif (display==7),			% keep in mind the last value of display 
  display=displayold;
 end

 if (linlogtfr==0),			% Lin/log scale of the tfr
  linlogstr='Change to a logarithmic scale';
 else
  linlogstr='Change to a linear scale';
 end;

 % Main menu
 choice=menu ('TFRQVIEW MENU :',...
	      'Change the display mode',...     
	      'Change the display layout',...   
	      'Change the color map',...        
	      'Change the sampling frequency',...
	      'Change the threshold',...
	      linlogstr,...
	      'Change the number of levels',...
	      'Grid',...        
	      'Save',...
	      'Save options',...	
	      'Print',...
	      'Close');

 if ((choice==9 | choice==11) & comp(1:2)~='PC'),
  if exist('workdir.mat'),
   load workdir;
  else  
   Path=pwd;
   str1=[' The current directory is ',Path];
   disp(str1); disp('');
   continue=0;
   while continue==0,
    continue=1;  
    str2=' Name of the directory (with full path) on which the files will be saved : ';
    disp(str2);
    PathWorkDir=input('       > ','s');
    if (exist(PathWorkDir)~=2),
     rep=input(' This directory doesn''t exist. Do you want to create it (y/n) ? ','s'); 
     if upper(rep)=='Y',
      eval(['!mkdir ',PathWorkDir]);
     else
      continue=0;
     end
    end
   end   
   if PathWorkDir(length(PathWorkDir))~='/',
    PathWorkDir=[PathWorkDir,'/'];
   end
   save workdir PathWorkDir;
  end
 end

 if (choice==1), 			% Change the display mode

  display=menu('DISPLAY MODE :',...
	       'Contour',...
	       'Imagesc',...
	       'Pcolor',...
	       'Surf',...
	       'Mesh',...
	       'Cancel');       

 elseif (choice==2), 			% Change the display layout

  layoutold=layout; indlo=1; layout=1;
  while layout~=6,			% While not close

   if ~iscolorbar,
    colorbarstr='Add a color bar';
   else
    colorbarstr='Remove the color bar';
   end
   if ~linlogspec,
    specstr='TFR + Spectrum (linear)';
    tfrsstr='TFR + Signal + Spectrum (linear)';
   else
    specstr='TFR + Spectrum (logarithmic)';
    tfrsstr='TFR + Signal + Spectrum (logarithmic)';
   end

   layout=menu('DISPLAY LAYOUT :',...
	       'Time-Frequency Representation',...
	       'TFR + Signal ',...
	       specstr,...
	       tfrsstr,...
	       colorbarstr,...
	       'Close');

   if (layout==3|layout==4),
    linlogspec=1-linlogspec;
   end
   if ((layout>=2 & layout<=4) & isempty(sig)),
    disp('Impossible action : the signal is unavailable');
   else
    access=2;
    linlog=linlogtfr+2*linlogspec;
    state=issig+2*iscolorbar;
    isgrid=isgridsig+2*isgridspe+4*isgridtfr;
    param = [display,linlog,threshold,levelnumb,Nf2,layout,access,state,fs,isgrid];
    map=colormap;
    if (nargin<=4),
     tfrview(tfr,sig,t,method,param,map);
    elseif (nargin==5),
     tfrview(tfr,sig,t,method,param,map,p1);
    elseif (nargin==6),
     tfrview(tfr,sig,t,method,param,map,p1,p2);
    elseif (nargin==7),
     tfrview(tfr,sig,t,method,param,map,p1,p2,p3);
    elseif (nargin==8),
     tfrview(tfr,sig,t,method,param,map,p1,p2,p3,p4);
    elseif (nargin==9),
     tfrview(tfr,sig,t,method,param,map,p1,p2,p3,p4,p5);
    end;

    axesh=get(gcf,'UserData');
    axcb=axesh(1); axsig=axesh(2); axspec=axesh(3); axtfr=axesh(4);

    if layout==1,			% Time-Frequency Representation
     axes(axsig);  cla reset; set(gca,'Visible','off');
     axes(axspec); cla reset; set(gca,'Visible','off');
     issig=0; isspec=0;
    elseif layout==2,			% TFR + Signal
     axes(axspec); cla reset; set(gca,'Visible','off');
     issig=1; isspec=0;
    elseif layout==3,			% TFR + spectrum
     axes(axsig);  cla reset; set(gca,'Visible','off');
     issig=0; isspec=1; 
    elseif layout==4,			% TFR + signal + spectrum
      issig=1; isspec=1; 
    elseif layout==5,			% Colorbar
     iscolorbar=1-iscolorbar;
    elseif ((layout==2|layout==3|layout==4) & isempty(sig)),
     disp('Unallowed action : the signal is anavailable.');
    end;

    lo(indlo)=layout;
    indlo=indlo+1;
   end
  end 
  while (layout==6|layout==5), 
   indlo=indlo-1; 
   if indlo>0,
    layout=lo(indlo); 
   else
    layout=layoutold;
   end
  end

 elseif (choice==3),			% Change the color map

  colmap=1;
  while colmap~=12,
   colmap=menu('COLOR MAP :',...
	       'hsv','jet','cool','bone','gray','hot','prism',...
	       'brighten','darken','permute','spin','Close');
   if colmap==1,
    colormap(hsv);
   elseif colmap==2,
    colormap(jet); 
   elseif colmap==3,
    colormap(cool);
   elseif colmap==4,
    colormap(bone);
   elseif colmap==5,
    colormap(gray);
   elseif colmap==6,
    colormap(hot);
   elseif colmap==7,
    colormap(prism);
   elseif colmap==8,
    brighten(0.25);
   elseif colmap==9,
    brighten(-0.25);
   elseif colmap==10,
    c = colormap; colormap(flipud(c));
   elseif colmap==11,
    spinmap;
   end
  end

 elseif (choice==4), 			% Change the sampling frequency

  fprintf(' Old sampling frequency: %f',fs); fsold=fs;
  fs=input(' New sampling frequency: ');
  if isempty(fs), fs=fsold; end

 elseif (choice==5),			% Change the threshold

  fprintf(' Old threshold: %f', threshold); throld=threshold;
  threshold=input(' New threshold: ');
  if isempty(threshold), threshold=throld; end

 elseif (choice==6),			% Change the lin/log scale of tfr

  linlogtfr=1-linlogtfr;

 elseif (choice==7),			% Change the number of levels

  fprintf(' Old number of levels: %f',levelnumb); levelold=levelnumb;
  levelnumb=input(' New number of levels: ');
  if isempty(levelnumb), levelnumb=levelold; end

 elseif (choice==8),			% Grids

  if ~issig & ~isspec,			% No signal and no spectrum
   isgridtfr=1-isgridtfr; axes(axtfr); grid

  elseif issig & ~isspec,		% A signal, no spectrum 
   Grid=1;
   while Grid~=3,
    if ~isgridsig,
     gridsigstr='Add a grid on the signal';
    else
     gridsigstr='Remove the grid on the signal';
    end
    if ~isgridtfr,
     gridtfrstr='Add a grid on the TFR';
    else
     gridtfrstr='Remove the grid on the TFR';
    end
    Grid=menu('GRID MENU :',gridsigstr,gridtfrstr,'Close');
    if Grid==1,
     isgridsig=1-isgridsig; axes(axsig); grid
    elseif Grid==2,
     isgridtfr=1-isgridtfr; axes(axtfr); grid  
    end
   end

  elseif ~issig & isspec,		% No signal, a spectrum
   Grid=1;
   while Grid~=3,
    if ~isgridspe,
     gridspestr='Add a grid on the spectrum';
    else
     gridspestr='Remove the grid on the spectrum';
    end
    if ~isgridtfr,
     gridtfrstr='Add a grid on the TFR';
    else
     gridtfrstr='Remove the grid on the TFR';
    end
    Grid=menu('GRID MENU :',gridspestr,gridtfrstr,'Close');
    if Grid==1,
     isgridspe=1-isgridspe; axes(axspec); grid
    elseif Grid==2,
     isgridtfr=1-isgridtfr; axes(axtfr); grid  
    end
   end

  else					% A signal and a spectrum
   Grid=1;
   while Grid~=4,
    if ~isgridsig,
     gridsigstr='Add a grid on the signal';
    else
     gridsigstr='Remove the grid on the signal';
    end
    if ~isgridspe,
     gridspestr='Add a grid on the spectrum';
    else
     gridspestr='Remove the grid on the spectrum';
    end
    if ~isgridtfr,
     gridtfrstr='Add a grid on the TFR';
    else
     gridtfrstr='Remove the grid on the TFR';
    end
    Grid=menu('GRID MENU :',gridsigstr,gridspestr,gridtfrstr,'Close');
    if Grid==1,
     isgridsig=1-isgridsig; axes(axsig); grid
    elseif Grid==2,
     isgridspe=1-isgridspe; axes(axspec); grid
    elseif Grid==3,
     isgridtfr=1-isgridtfr; axes(axtfr); grid
    end
   end
  end

 elseif (choice==9),			% Save the parameters

  f=freq*fs;
  Nmethod=length(method);
  if comp(1:2)=='PC',
   namedflt=[method(4:Nmethod),num2str(Nsig),'.mat'];
   [name,PathWorkDir] = uiputfile(namedflt, 'Save As');
  else
   namedflt=[method(4:Nmethod),num2str(Nsig)];
   nameStr=[' Name of the MAT file [',namedflt,'] : '];
   name=input(nameStr,'s'); 
   while (length(name)>8),
    disp(' The name must have less than 8 characters');
    name=input(nameStr,'s'); 
   end
   if isempty(name),
    name=namedflt;
   end
  end
  access=0; 
  linlog=linlogtfr+2*linlogspec;
  state=issig+2*iscolorbar;
  isgrid=isgridsig+2*isgridspe+4*isgridtfr;
  param = [display,linlog,threshold,levelnumb,Nf2,layout,access,state,fs,isgrid];
  map=colormap;
  if (nargin<=4),
   TfrQView=['tfrqview(tfr,sig,t,method)'];
   TfrView =['clf;tfrview(tfr,sig,t,method,param,map)'];
   eval(['save ',PathWorkDir,name,...
' tfr sig t f fs method param map TfrView TfrQView']);
  elseif (nargin==5),
   TfrQView=['tfrqview(tfr,sig,t,method,p1)'];
   TfrView =['clf; tfrview(tfr,sig,t,method,param,map,p1)'];
   eval(['save ',PathWorkDir,name,...
' tfr sig t f fs method param map p1 TfrView TfrQView']);
  elseif (nargin==6),
   TfrQView=['tfrqview(tfr,sig,t,method,p1,p2)'];
   TfrView =['clf; tfrview(tfr,sig,t,method,param,map,p1,p2)'];
   eval(['save ',PathWorkDir,name,...
' tfr sig t f fs method param map p1 p2 TfrView TfrQView']);
  elseif (nargin==7),
   TfrQView=['tfrqview(tfr,sig,t,method,p1,p2,p3)'];
   TfrView =['clf; tfrview(tfr,sig,t,method,param,map,p1,p2,p3)'];
   eval(['save ',PathWorkDir,name,...
' tfr sig t f fs method param map p1 p2 p3 TfrView TfrQView']);
  elseif (nargin==8),
   TfrQView=['tfrqview(tfr,sig,t,method,p1,p2,p3,p4)'];
   TfrView =['clf; tfrview(tfr,sig,t,method,param,map,p1,p2,p3,p4)'];
   eval(['save ',PathWorkDir,name,...
' tfr sig t f fs method param map p1 p2 p3 p4 TfrView TfrQView']);
  elseif (nargin==9),
   TfrQView=['tfrqview(tfr,sig,t,method,p1,p2,p3,p4,p5)'];
   TfrView =['clf; tfrview(tfr,sig,t,method,param,map,p1,p2,p3,p4,p5)'];
   eval(['save ',PathWorkDir,name,...
' tfr sig t f fs method param map p1 p2 p3 p4 p5 TfrView TfrQView']);
  end;
  disp(' ');
  SaveTxt1=[' The file is saved in the directory ',PathWorkDir];
  SaveTxt2=['  under the name ',name,'.mat'];
  disp(SaveTxt1); disp(SaveTxt2);

  disp(' If you want to find again the exact layout of this screen, do');
  Txt1=['	load ',PathWorkDir,name,'; eval(TfrView);'];
  disp(Txt1);   
  disp(' If you want to restart the display session under tfrqview, do');
  Txt2=['	load ',PathWorkDir,name,'; eval(TfrQView);'];
  disp(Txt2);   

 elseif (choice==10),			% Save options

  save options fs threshold linlogtfr linlogspec levelnumb colmap ...
       display layout iscolorbar isgridsig ...
       isgridspe isgridtfr issig isspec;
  disp(' '); disp(' Options saved');
 
 elseif (choice==11),			% Print the current figure

  Nmethod=length(method);
  if comp(1:2)=='PC',
   namedflt=[method(4:Nmethod),num2str(Nsig),'.eps'];
   [name,PathWorkDir] = uiputfile(namedflt, 'Save As');
  else
   namedflt=[method(4:Nmethod),num2str(Nsig)];
   nameStr=[' Name of the EPS file [',namedflt,'] : '];
   name=input(nameStr,'s'); 
   while (length(name)>8),
    disp('The name must have less than 8 characters');
    name=input(nameStr,'s'); 
   end
   if isempty(name),
    name=namedflt;
   end
  end
  eval(['print -deps ',PathWorkDir,name]); disp(' ');
  PrintTxt1=[' The file is saved in the directory ',PathWorkDir];
  if comp(1:2)=='PC',
   PrintTxt2=['  under the name ',name];
  else
   PrintTxt2=['  under the name ',name,'.eps'];
  end 
  disp(PrintTxt1); disp(PrintTxt2);

 end;

 if display~=7,
  displayold=display;
 end

end

if (comp(1:2)~='PC' & exist('workdir.mat')),
 !rm workdir.mat
end

disp(' ');
