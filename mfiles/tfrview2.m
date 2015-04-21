function tfrview(tfr,sig,t,method,param,map,p1,p2,p3,p4,p5);
%TFRVIEW Visualization of time-frequency representations.
%	TFRVIEW(TFR,SIG,T,METHOD,PARAM,MAP,P1,P2,P3,P4,P5) 
%	allows to visualize a time-frequency representation.
%	TFRVIEW is called through TFRQVIEW from any TFR* function.
%
%	TFR    : time-frequency representation.
%	SIG    : signal in the time-domain.  
%	T      : time instants.
%	METHOD : chosen representation (name of the corresponding M-file)
%	PARAM  : visualization parameter vector :
%	 PARAM = [DISPLAY LINLOG THRESHOLD LEVNUMB NF2 LAYOUT 
%				ACCESS STATE FS ISGRID] where
%	  DISPLAY=1..5 for contour, imagesc, pcolor, surf or mesh
%	  LINLOG=0/1 for linearly/logarithmically spaced levels
%	  THRESHOLD  is the visualization threshold, in % 
%	  LEVELNUMB  is the number of levels used with contour
%	  NF2        is the number of frequency bins displayed
%	  LAYOUT determines the layout of the figure : TFR alone (1),  
%		TFR and SIG (2), TFR and spectrum (3), TFR and SIG and 
%		spectrum (4), add/remove the colorbar (5)
%	  ACCESS depends on the way you access to tfrview : from the 
%		command line (0) ; from tfrqview, except after a
%		change in the sampling frequency or in the layout (1) ; 
%		from tfrqview, after a change in the layout (2) ; 
%		from tfrqview, after a change in the sampling frequency (3)
%	  STATE depends on the signal/colorbar presence : 
%		no signal, no colorbar (0) ; signal, no colorbar (1) ;
%		no signal, colorbar (2) ; signal and colorbar (3)
%	  FS is the sampling frequency (may be set to 1.0)
%	  ISGRID depends on the grids' presence : 
%			isgrid=isgridsig+2*isgridspe+4*isgridtfr
%		where isgridsig=1 if a grid is present on the signal
%		and =0 if not, and so on
%	MAP    :  selected colormap. 
%	P1..P5 : parameters of the representation. Run the file 
%		 TFRPARAM(METHOD) to know the meaning of P1..P5.  		 
%
%	TFRVIEW is called through TFRQVIEW by any file TFR*.
%
%			Use TFRQVIEW preferably.
%			------------------------
%
%	See also TFRQVIEW, TFRPARAM.

%	F. Auger, July 1994, July 1995 - 
%	O. Lemoine, October-November 1995, May-June 1996. 
%	Copyright (c) CNRS - France 1996. 
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

if ( nargin < 6 ),
 error ('at least 6 parameters are required');
end;
[tfrrow,tfrcol] = size(tfr);
[trow,tcol] = size(t);
[Nsig,Ncol]=size(sig);
maxi=max(max(tfr)); 

% Extraction of the elements of param
display   = param(1); 
linlog    = param(2); 
threshold = param(3); 
levelnumb = param(4);
Nf2       = param(5);
layout	  = param(6);
access    = param(7);
state     = param(8);
fs        = param(9);
isgrid    = param(10);

if fs<1000,
 unitHz=1;
elseif (fs>=1e3 & fs<1e6),
 fs=fs/1e3;
 unitHz=2;
elseif (fs>=1e6),
 fs=fs/1e6;
 unitHz=3;
end

linlogtfr=rem(linlog,2);
linlogspec=(linlog-linlogtfr)/2;
issig=rem(state,2);
iscolorbar=(state-issig)/2;
isgridsig=rem(isgrid,2);
tempo=(isgrid-isgridsig)/2;
isgridspe=rem(tempo,2);
isgridtfr=(tempo-isgridspe)/2;


% Computation of isaffine and freq (vector of frequency samples) 
if strcmp(method,'TFRASPW' ) | strcmp(method,'TFRSCALO') | ...
   strcmp(method,'TFRDFLA' ) | strcmp(method,'TFRSPAW' ) | ...
   strcmp(method,'TFRUNTER') | strcmp(method,'TFRBERT' ),
 freq=eval(['p',num2str(nargin-6)]);   	% last input argument is freqs.
 isaffine=1;
 if display==2,				% imagesc do not allow
  display=3;				% non linear scales for the axes.
  disp('Imagesc does not support non-linear scales for axes. We use pcolor instead');
 end
else
 isaffine=0;
 freq=(0.5*(0:Nf2-1)/Nf2);
end
freqr=freq*fs;
ts=t/fs;


% Update variables mini, levels, Strlinlog according to linlog 
if ~linlogtfr,
 if (display==4|display==5),
  mini=min(min(tfr));
 else
  mini=max(min(min(tfr)),maxi*threshold/100.0);
 end
 levels=linspace(mini,maxi,levelnumb+1); 
 Strlinlog=', lin. scale';
else
 mini=max(min(min(tfr)),maxi*threshold/100.0);
 levels=logspace(log10(mini),log10(maxi),levelnumb+1);
 Strlinlog=', log. scale';
end;


% Initialization of the axes
if (access==1|access==2|access==3),			% From the main menu of tfrqview
 axesh  = get(gcf,'UserData');
 axcb   = axesh(1);
 axsig  = axesh(2);
 axspec = axesh(3);
 axtfr  = axesh(4);
elseif (access==0),		% From the save option of tfrqview
 axcb   = axes('Units','normal','Visible','off','Box','On');
 axsig  = axes('Units','normal','Visible','off','Box','On');
 axspec = axes('Units','normal','Visible','off','Box','On');
 axtfr  = axes('Units','normal','Visible','off','Box','On');
end

if (access==0|access==2|access==3),

 % Test of analycity and computation of spec
 if ~isempty(sig),
  spec=zeros(2*Nf2,Ncol);
  for k=1:Ncol,
   isana=1; alpha=2; Lt=max(t)-min(t)+1;   	
   if 2*Nf2>=Lt,
    spec(:,k)=abs(fft(sig(min(t):max(t),k),2*Nf2)).^2; 
   else
    sp=abs(fft(sig(min(t):max(t),k))).^2; 
    fr1=(0.5*(0:Lt-1)/Lt)*fs;
    fr2=(0.5*(0:2*Nf2-1)/2/Nf2)*fs;
    spec(:,k)=interp1(fr1,sp,fr2)';
   end
   spec1=sum(spec(1:Nf2,k));
   spec2=sum(spec(Nf2+1:2*Nf2,k));
   if spec2>spec1/10,
    isana=0;
    if ~isreal(sig(min(t):max(t),k)),
     alpha=1;
    end
   end
  end
 end

 if layout==1,			% Time-Frequency Representation

  set(axtfr,'Position',[0.13 0.11 0.775 0.815]);
  issig=0; isspec=0;

 elseif layout==2,			% TFR + Signal

  set(axtfr,'Position',[.1 .1 .83 .55]);
  set(axsig,'Position',[.1 .72 .83 .2]); 
  axes(axsig);  
  plot(ts,real(sig(t,:)));
  set(gca,'xticklabel',[]);
  if comp(1:2)=='PC', set(gca,'fontsize',10); end
  ylabel('Real part');
  title('Signal in time');
  Min=min(min(real(sig))); Max=max(max(real(sig)));
  axis([min(ts) max(ts) Min Max]);
  issig=1; isspec=0;

 elseif layout==3,			% TFR + spectrum

  axes(axspec); 
  set(gca,'Position',[.1 .1 .18 .8]);
  if isaffine, 
   f1=freqr(1); f2=freqr(Nf2); df=f2-f1; 
   Nf4=round((Nf2-1)*fs/(2*df))+1;
   for k=1:Ncol,
    spec(1:alpha*Nf4,k)=abs(fft(sig(min(t):max(t),k),alpha*Nf4)).^2; 
   end
   spec=spec((round(f1*2*(Nf4-1)/fs)+1):(round(f1*2*(Nf4-1)/fs)+Nf2),:);
   freqs=linspace(f1,f2,Nf2);
  else
   freqs=freqr;
   spec=spec(1:Nf2,:);
  end
  Maxsp=max(max(spec));
  if linlogspec,
   plot(freqs,spec);
   if comp(1:2)=='PC', set(gca,'fontsize',10); end
   title('Linear scale'); 
   set(gca,'ytick',[0 round(Maxsp/2) round(Maxsp)]);
  else
   semilogy(freqs,spec);
   if comp(1:2)=='PC', set(gca,'fontsize',10); end
   title('Log. scale [dB]');
   set(gca,'ytick',[fix(10^(log10(Maxsp)/2));fix(Maxsp)]);
   str1=['   ']; str2=['   '];
   st1=num2str(fix(10*log10(Maxsp)/2)); str1(1:length(st1))=st1;
   st2=num2str(fix(10*log10(Maxsp)));   str2(1:length(st2))=st2; 
   set(axspec,'yticklabel',[str1;str2]);
  end
  xlabel('Energy spectral density');
  Nsp=length(spec); 
  set(gca,'Xlim',[freqs(1) freqs(Nsp)]);
  set(gca,'Ylim',[0 Maxsp*1.2]);
  set(gca,'xticklabel',[],'view',[-90 90]);

  set(axtfr,'Position',[.36 .1 .57 .8]);
  issig=0; isspec=1; 

 elseif layout==4,			% TFR + signal + spectrum
 
  set(axsig,'Position',[.33 .73 .6 .2]);
  axes(axsig); 
  plot(ts,real(sig(t,:)));
  if comp(1:2)=='PC', set(gca,'fontsize',10); end
  ylabel('Real part');
  title('Signal in time');
  set(gca,'xticklabel',[]);
  Min=min(min(real(sig))); Max=max(max(real(sig)));
  axis([min(ts) max(ts) Min Max]);

  set(axspec,'Position',[.1 .1 .15 .55]);
  axes(axspec); 
  if isaffine, 
   f1=freqr(1); f2=freqr(Nf2); df=f2-f1; 
   Nf4=round((Nf2-1)*fs/(2*df))+1;
   for k=1:Ncol,
    spec(1:alpha*Nf4,k)=abs(fft(sig(min(t):max(t),k),alpha*Nf4)).^2; 
   end
   spec=spec((round(f1*2*(Nf4-1)/fs)+1):(round(f1*2*(Nf4-1)/fs)+Nf2),:);
   freqs=linspace(f1,f2,Nf2);
  else
   freqs=freqr;
   spec=spec(1:Nf2,:);
  end
  Maxsp=max(max(spec));
  if linlogspec,
   plot(freqs,spec);
   if comp(1:2)=='PC', set(gca,'fontsize',10); end
   title('Linear scale');  
   set(axspec,'ytick',[0 round(Maxsp/2) round(Maxsp)]);
  else
   semilogy(freqs,spec); 
   if comp(1:2)=='PC', set(gca,'fontsize',10); end
   title('Log. scale [dB]');
   set(gca,'ytick',[fix(10^(log10(Maxsp)/2));fix(Maxsp)]);
   str1=['   ']; str2=['   '];
   st1=num2str(fix(10*log10(Maxsp)/2)); str1(1:length(st1))=st1;
   st2=num2str(fix(10*log10(Maxsp)));   str2(1:length(st2))=st2; 
   set(axspec,'yticklabel',[str1;str2]);
  end
  xlabel('Energy spectral density');
  Nsp=length(spec); 
  set(gca,'Xlim',[freqs(1) freqs(Nsp)]);
  set(gca,'Ylim',[0 Maxsp*1.2]);
  set(gca,'xticklabel',[],'view',[-90 90]);

  set(axtfr,'Position',[.33 .1 .6 .55]);
  issig=1; isspec=1;

 elseif layout==5,			% Colorbar

  if ~iscolorbar,                     % there is no color bar already 
   pos=get(axtfr,'Position');
   set(axtfr,'Position',[pos(1) pos(2) pos(3)-0.03 pos(4)]);
   set(axcb,'Position',[pos(1)+pos(3)-0.01 pos(2) 0.01 pos(4)]);
   axes(axcb); 
   Ncolors=length(colormap); 
   [cmin,cmax]=caxis; colorvect=linspace(cmax,cmin,Ncolors);
   imagesc(colorvect'); axis('off'); 
   if issig,                          % there is a signal
    possig=get(axsig,'Position');
    set(axsig,'Position',[possig(1) possig(2) possig(3)-0.03 possig(4)]);
   end 
  else                                % there is already a color bar
   pos=get(axtfr,'Position');
   set(axtfr,'Position',[pos(1) pos(2) pos(3)+0.03 pos(4)]);
   axes(axcb); cla reset; set(gca,'Visible','off'); 
   if issig,                          % there is a signal
    possig=get(axsig,'Position');
    set(axsig,'Position',[possig(1) possig(2) possig(3)+0.03 possig(4)]);
   end
  end
  iscolorbar=1-iscolorbar;
  isspec=[];
 end;

 if (layout<5 & iscolorbar),		% Updating of the colorbar
  pos=get(axtfr,'Position');
  set(axtfr,'Position',[pos(1) pos(2) pos(3)-0.03 pos(4)]);
  set(axcb,'Position',[pos(1)+pos(3)-0.01 pos(2) 0.01 pos(4)]);
  axes(axcb); 
  Ncolors=length(colormap); 
  [cmin,cmax]=caxis; colorvect=linspace(cmax,cmin,Ncolors);
  imagesc(colorvect'); axis('off'); 
  if issig,                          % there is a signal
   possig=get(axsig,'Position');
   set(axsig,'Position',[possig(1) possig(2) possig(3)-0.03 possig(4)]);
  end 
 end

 if (layout>5),	                  	% Avoid warning messages occuring with Matlab 5
  isspec=[];
 end

 set(gcf,'UserData',[axcb,axsig,axspec,axtfr]);
 axes(axtfr);
end

if (access==0|access==1|access==3),

 % Display the tfr
 if (tcol==1),
   plot(freqr,tfr(1:Nf2));
 else
   if strcmp(computer,'MAC'),
     tfr=flipud(tfr(1:Nf2,:));
   else
     tfr=tfr(1:Nf2,:);
   end;
 
   indmin=find(tfr<mini);
   tfr(indmin)=mini*ones(1,length(indmin));
 
   indmax=find(tfr>maxi);
   tfr(indmax)=maxi*ones(1,length(indmax));
   
   if display==1,
    vnum=version;
    if vnum(1)=='5'
     contour(ts,freqr,tfr,levels);
    else
     contour(tfr,levels,ts,freqr);
    end;
   elseif display==2,
    if linlogtfr==0,
     imagesc(ts,freqr,tfr);axis('xy');
    else
     imagesc(ts,freqr,log10(tfr));axis('xy');
    end
   elseif display==3,
    if linlogtfr==0,
     pcolor(ts,freqr,tfr); shading interp;
    else
     pcolor(ts,freqr,log10(tfr)); shading interp;
    end
   elseif display==4,
    if linlogtfr==0,
     surf(ts,freqr,tfr); shading interp;
     if comp(1:2)=='PC', set(gca,'fontsize',10); end
     zlabel('Amplitude');
     axis([ts(1) ts(tcol) freqr(1) freqr(Nf2) mini maxi]);
    else
     surf(ts,freqr,log10(tfr)); shading interp;
     if comp(1:2)=='PC', set(gca,'fontsize',10); end
     zlabel('Positive values');
     axis([ts(1) ts(tcol) freqr(1) freqr(Nf2) log10(mini) log10(maxi)]);
    end
   elseif display==5,
    if linlogtfr==0,
     mesh(ts,freqr,tfr); shading interp;
     if comp(1:2)=='PC', set(gca,'fontsize',10); end
     zlabel('Amplitude');
     axis([ts(1) ts(tcol) freqr(1) freqr(Nf2) mini maxi]);
    else
     mesh(ts,freqr,log10(tfr));shading interp;
     if comp(1:2)=='PC', set(axtfr,'fontsize',10); end
     zlabel('Positive values');
     axis([ts(1) ts(tcol) freqr(1) freqr(Nf2) log10(mini) log10(maxi)]);
    end
   end
 end;
 
 % Define the title and check the input arguments depending on 'method'

 if comp(1:2)=='PC',
  set(axsig  ,'fontsize',10); 
  set(axspec ,'fontsize',10); 
  set(axtfr  ,'fontsize',10); 
 end

 method=method(4:length(method));

 if nargin==6,
  title([method, Strlinlog,...
         ', Threshold=',num2str(threshold),'%']);

 elseif strcmp(method,'WV') | strcmp(method,'MH') | ...
   strcmp(method,'PAGE'),
  title([method,', Nf=',num2str(Nf2), Strlinlog,...
        ', Threshold=',num2str(threshold),'%']);

 elseif strcmp(method,'PWV'  )|strcmp(method,'PMH'  )| ...
        strcmp(method,'SP'   )|strcmp(method,'PPAGE')| ...
        strcmp(method,'RSP'  )|strcmp(method,'RPPAG')| ...
        strcmp(method,'RPWV' )|strcmp(method,'RPMH' ),
  h=p1;[hrow,hcol]=size(h); Lh=(hrow-1)/2;
  if (hcol~=1)|(rem(hrow,2)==0),
   error('h must be a smoothing window with odd length'); end;
  title([method, ', Lh=',num2str(Lh), ', Nf=',num2str(Nf2),...
        Strlinlog,', Threshold=',num2str(threshold),'%']);
 
 elseif strcmp(method,'STFT'),
  h=p1;[hrow,hcol]=size(h); Lh=(hrow-1)/2;
  if (hcol~=1)|(rem(hrow,2)==0),
   error('h must be a smoothing window with odd length'); end;
  title(['|',method, '|^2, Lh=',num2str(Lh),...
        ', Nf=',num2str(Nf2), Strlinlog,', Thld=',...
        num2str(threshold),'%']);
 
 elseif strcmp(method,'SPWV' ) | strcmp(method,'MHS'  )| ...
        strcmp(method,'RSPWV') | strcmp(method,'RMHS' )| ...
        strcmp(method,'ZAM'  ) | strcmp(method,'RIDBN')|...
        strcmp(method,'BJ'   ) | strcmp(method,'RIDB' )| ...
        strcmp(method,'RIDH' ) | strcmp(method,'RIDT' ),
  g=p1; [grow,gcol]=size(g); Lg=(grow-1)/2;
  if (gcol~=1)|(rem(grow,2)==0),
   error('g must be a smoothing window with odd length'); end;
  h=p2; [hrow,hcol]=size(h); Lh=(hrow-1)/2; 
  if (hcol~=1)|(rem(hrow,2)==0),
   error('h must be a smoothing window with odd length'); end;
  title([method,', Lg=',num2str(Lg),', Lh=',num2str(Lh),...
         ', Nf=',num2str(Nf2),Strlinlog,...
         ', Threshold=',num2str(threshold),'%']);
 
 elseif strcmp(method,'MMCE'),
  h=p1;[hrow,hcol]=size(h); Lh=(hrow-1)/2;
  if (rem(hrow,2)==0),
   error('h must be a smoothing window with odd length'); end;
  title([method, ', Lh=',num2str(Lh), ', Nf=',num2str(Nf2),...
        Strlinlog,', Threshold=',num2str(threshold),'%']);
 
 elseif strcmp(method,'CW' ) | strcmp(method,'BUD'),
  g=p1; [grow,gcol]=size(g); Lg=(grow-1)/2;
  if (gcol~=1)|(rem(grow,2)==0),
   error('g must be a smoothing window with odd length'); end;
  h=p2; [hrow,hcol]=size(h); Lh=(hrow-1)/2; 
  if (hcol~=1)|(rem(hrow,2)==0),
   error('h must be a smoothing window with odd length'); end;
  sigma=p3;
  title([method,', Lg=',num2str(Lg),', Lh=',num2str(Lh),...
         ' sigma=',num2str(sigma),', Nf=',num2str(Nf2),...
         Strlinlog, ', Threshold=',num2str(threshold),'%']);
 
 elseif strcmp(method,'GRD')
  g=p1; [grow,gcol]=size(g); Lg=(grow-1)/2;
  if (gcol~=1)|(rem(grow,2)==0),
   error('g must be a smoothing window with odd length'); end;
  h=p2; [hrow,hcol]=size(h); Lh=(hrow-1)/2; 
  if (hcol~=1)|(rem(hrow,2)==0),
   error('h must be a smoothing window with odd length'); end;
  title([method,', Lg=',num2str(Lg),', Lh=',num2str(Lh),...
         ', rs =',num2str(p3), ', M/N =',num2str(p4), ...
         ', Nf =',num2str(Nf2), ...
         Strlinlog, ', Threshold=',num2str(threshold),'%']);
 
 elseif strcmp(method,'MSC' ) | strcmp(method,'RMSC' )
  f0T=p1; if (f0T<=0), error('f0T must be positive'); end;
  title([method, ', f0T=',num2str(f0T), ', Nf=',num2str(Nf2),...
        Strlinlog, ', Threshold=',num2str(threshold),'%']);
 
 elseif strcmp(method,'RGAB' )
  Nh=p1; if (Nh<=0), error('Nh must be positive'); end;
  title([method, ', Nh=',num2str(Nh), ', Nf=',num2str(Nf2),...
        Strlinlog, ', Threshold=',num2str(threshold),'%']);
 
 
 elseif strcmp(method,'DFLA' ) | strcmp(method,'UNTER' )| ...
        strcmp(method,'BERT' ),
  N=p1; 
  if (N<=0),                error('N must be positive'); end;
  title([method, ', N=',num2str(N), Strlinlog, ', Threshold=',...
         num2str(threshold), '%']);  
 
 elseif strcmp(method,'SCALO'),
  Nh0=p1; N=p2; 
  if (Nh0<0),               error('Nh0 must be positive'); end;
  if (N<=0),                error('N must be positive'); end;
  if (Nh0>0), 
   title([method, ', Morlet wavelet, Nh0=', num2str(Nh0), ...
         ', N=',num2str(N), Strlinlog, ', Thld=',...
         num2str(threshold), '%']);  
  else 
   title([method, ', Mexican hat, N=',num2str(N), Strlinlog, ', Thld=',...
         num2str(threshold), '%']);  
  end
 
 elseif strcmp(method,'ASPW'),
  Nh0=p1; Ng0=p2; N=p3; 
  if (Nh0<0),               error('Nh0 must be positive'); end;
  if (Ng0<0),               error('Ng0 must be positive'); end;
  if (N<=0),                error('N must be positive'); end;
  if (Nh0>0), 
   title([method, ', Morlet wlt, Nh0=', num2str(Nh0), ', Ng0=',...
         num2str(Ng0), ', N=',num2str(N), Strlinlog, ', Thld=',...
         num2str(threshold), '%']);  
  else 
   title([method, ', Mexican hat, Ng0=',num2str(Ng0),...  
         ', N=',num2str(N), Strlinlog, ', Thld=',... 
         num2str(threshold), '%']);  
  end
 
 elseif strcmp(method,'SPAW'),
  K=p1; Nh0=p2; Ng0=p3; N=p4; 
  if (Nh0<0),               error('Nh0 must be positive'); end;
  if (Ng0<0),               error('Ng0 must be positive'); end;
  if (N<=0),                error('N must be positive'); end;
  if (Nh0>0), 
   title([method, ', K=', num2str(K), ', Morlet wlt, Nh0=',...
         num2str(Nh0), ', Ng0=',...
         num2str(Ng0), ', N=',num2str(N), Strlinlog, ', Thld=',...
         num2str(threshold), '%']);  
  else 
   title([method, ', K=', num2str(K), ', Mexican hat, Ng0=',...
         num2str(Ng0),', N=',num2str(N), Strlinlog, ', Thld=',...
         num2str(threshold), '%']);  
  end
 
 
 elseif strcmp(method,'GABOR'),
  N=p1; Q=p2; h=p3; [hrow,hcol]=size(h); Lh=(hrow-1)/2;
  if (hcol~=1)|(rem(hrow,2)==0),
   error('h must be a smoothing window with odd length'); end;
  title([method, ', Lh=',num2str(Lh), ', Nf=',...
         num2str(Nf2),', N=',num2str(N),', Q=',num2str(Q),...
         Strlinlog,', Thld=',num2str(threshold),'%']);
 
 end;
end

if unitHz==1,
 xlabel('Time [s]');
 ylabel('Frequency [Hz]');
elseif unitHz==2,
 xlabel('Time [ms]');
 ylabel('Frequency [kHz]');
elseif unitHz==3,
 xlabel('Time [µs]');
 ylabel('Frequency [MHz]');
end

if (access==1),                 % Avoid warning message occuring Matlab 5
 isspec=[];
end;

if (isgridsig & issig),		% Updating of the grids
 axes(axsig); grid on	
elseif (~isgridsig & issig),
 axes(axsig); grid off
end 
if (isgridspe & isspec),
 axes(axspec); grid on
elseif (~isgridspe & isspec),
 axes(axspec); grid off
end 

if access~=2,
 if (isgridtfr),			% upating of the grid on the tfr
  axes(axtfr); grid on  
 elseif (~isgridtfr),
  axes(axtfr); grid off
 end 
end

if (access==0),
 colormap(map);
end
