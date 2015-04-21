function tfrview(tfr,sig,t,method,param,p1,p2,p3,p4,p5);
%TFRVIEW Visualization of time-frequency representations.
%        TFRVIEW(TFR,SIG,T,METHOD,PARAM,P1,P2,P3,P4,P5) 
%        allows to visualize a time-frequency representation.
%        TFRVIEW is called through TFRQVIEW from any TFR* function.
%
%        TFR    : time-frequency representation.
%        SIG    : signal in the time-domain.  
%        T      : time instants.
%        METHOD : chosen representation (name of the corresponding M-file)
%        PARAM  : visualization parameter vector :
%         PARAM  = [DISPLAY LINLOG THRESHOLD LEVNUMB NF2 LAYOUT FS ISGRID] where
%         DISPLAY=1..5 for contour, imagesc, pcolor, surf or mesh
%         LINLOG =0/1 for linearly/logarithmically spaced levels
%                  
%         THRESHOLD  is the visualization threshold, in % 
%         LEVELNUMB  is the number of levels used with contour
%         NF2        is the number of frequency bins displayed
%         LAYOUT determines the layout of the figure : TFR alone (1),  
%                TFR and SIG (2), TFR and spectrum (3), TFR and SIG and 
%                spectrum (4), add 4 if you want a colorbar
%         FS     is the sampling frequency (may be set to 1.0)
%         ISGRID depends on the grids' presence : 
%                isgrid=isgridsig+2*isgridspec+4*isgridtfr
%                where isgridsig=1 if a grid is present on the signal
%                and =0 if not, and so on
%         fmin   smallest normalized frequency
%         fmax   highest  normalized frequency
%        P1..P5: parameters of the representation. Run the file 
%                TFRPARAM(METHOD) to know the meaning of P1..P5.  		 
%
%	TFRVIEW is called through TFRQVIEW by any file TFR*.
%
%			Use TFRQVIEW preferably.
%			------------------------
%
%	See also TFRQVIEW, TFRPARAM.

%	F. Auger, July 1994, July 1995 - 
%	O. Lemoine, October-November 1995, May-June 1996. 
%       F. Auger, May 1998.
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

comp=computer; % so as to know the running computer
MatlabVersion=version; MatlabVersion=str2num(MatlabVersion(1));
% I hope that future Matlab versions will be more compatible
if (MatlabVersion==4),
 TickLabelStr='TickLabels';
elseif (MatlabVersion>=5),
 TickLabelStr='Ticklabel';
else
 error('unsupported matlab version. please send an email.'); 
end;
 
if ( nargin < 5 ),
 error ('at least 5 parameters are required');
end;

[tfrrow,tfrcol] = size(tfr); % the size of tfr
[trow,tcol]     = size(t);   % the size of t
[Nsig,Ncol]     = size(sig); % the size of sig
maxi=max(max(tfr)); 

% Extraction of the elements of param
display     = param( 1);      % contour, imagesc, pcolor, surf, or mesh
linlog      = param( 2);      % linear or logarithmic scale
threshold   = param( 3);      % visualization threshold
levelnumb   = param( 4);      % number of levels
Nf2         = param( 5);      % number of frequency points
layout      = param( 6);      % figure layout
fs          = param( 7);      % sampling frequency
isgrid      = param( 8);      % grid(s)
fmin        = param( 9);      % smallest displayed frequency 
fmax        = param(10);      % highest displayed frequency

if (fmin>=fmax),
 error('fmin should be lower than fmax');
elseif (fmax>0.5),
 error('fmax is a normalized frequency, and should be lower than 0.5');
end

if fs<1000, unitHz=1;                            % display in  s and  Hz
elseif (fs>=1e3 & fs<1e6), fs=fs/1e3; unitHz=2;  % display in ms and kHz
elseif (fs>=1e6), fs=fs/1e6; unitHz=3;           % display in us and MHz
end

linlogtfr=rem(linlog,2);
linlogspec=rem(linlog-linlogtfr,4)/2;
sigenveloppe=rem(linlog-linlogtfr-linlogspec*2,8)/4;

issig=rem(layout-1,2);
isspec=rem(layout-1-issig,4)/2;
iscolorbar=rem(layout-1-issig-isspec*2,8)/4;
if isempty(sig),                        % I can't make miracles
 issig=0; isspec=0; layout=issig+isspec*2+1;
else  
 layout=layout-4*iscolorbar;
end;

isgridsig =rem(isgrid,2);
isgridspec=rem(isgrid-isgridsig,4)/2;
isgridtfr =rem(isgrid-isgridsig-isgridspec*2,8)/4;

% Computation of isaffine and freq (vector of frequency samples) 
if istfraff(method),
 freq=eval(['p',num2str(nargin-5)]);   	% last input argument is freqs.
 isaffine=1;
 if display==2,	                        % imagesc not allowed
  display=3;                            % non linear scales for the axes.
  disp('Imagesc does not support non-linear scales for axes. We use pcolor instead');
 end
else
 isaffine=0;                            % weyl-heisenberg group of distributions
 freq=(0.5*(0:Nf2-1)/Nf2);              % dispaly only the positive frequencies
end
freqr=freq*fs; ts=t/fs;                 % real time and frequency


% Update variables mini, levels, LinLogStr according to linlog 
if ~linlogtfr,
 if (display==4)|(display==5),          % surf or mesh
  mini=min(min(tfr));
 else
  mini=max(min(min(tfr)),maxi*threshold/100.0);
 end
 levels=linspace(mini,maxi,levelnumb+1); 
 LinLogStr=', lin. scale';
else
 mini=max(min(min(tfr)),maxi*threshold/100.0);
 levels=logspace(log10(mini),log10(maxi),levelnumb+1);
 LinLogStr=', log. scale';
end;

% Initialization of the current figure
zoom off; clf; 
set(gcf,'Resize','On','NextPlot','Add');

% Initialization of the axes
if iscolorbar,
 axcb   = axes('Units','normal','Visible','off','Box','On');
 set(gcf,'UserData',[get(gcf,'UserData') axcb]);
end;

if issig,
 axsig  = axes('Units','normal','Visible','off','Box','On');
 if comp(1:2)=='PC', set(axsig ,'fontsize',10); end
 set(gcf,'UserData',[get(gcf,'UserData') axsig]);
end;

if isspec,
 axspec = axes('Units','normal','Visible','off','Box','On');
 if comp(1:2)=='PC', set(axspec,'fontsize',10); end;
 set(gcf,'UserData',[get(gcf,'UserData') axspec]);
end;

axtfr  = axes('Units','normal','Visible','off','Box','On');
if comp(1:2)=='PC', set(axtfr ,'fontsize',10); end
set(gcf,'UserData',[get(gcf,'UserData') axtfr]);

 % Test of analycity and computation of spec
 if ~isempty(sig),
  for k=1:Ncol,
   isana=1; alpha=2; Lt=max(t)-min(t)+1;   	
   if 2*Nf2>=Lt,
    spec(:,k)=abs(fft(sig(min(t):max(t),k),2*Nf2)).^2; 
   else
    % modifications :  F. Auger (fog), 30/11/97
    Nb_tranches_fog = floor(Lt/(2*Nf2));
    % fprintf('%f \n',Nb_tranches_fog);
    spec(:,k)=zeros(2*Nf2,1);
    for Num_tranche_fog=0:Nb_tranches_fog-1,
     % fprintf('%f \n',Num_tranche_fog);
     spec(:,k)=spec(:,k)+abs(fft(sig(min(t)+2*Nf2*Num_tranche_fog+(0:2*Nf2-1),k))).^2;
    end;

    if (Lt>Nb_tranches_fog*2*Nf2),
     spectre_fog=fft(sig(min(t)+tfrrow*Nb_tranches_fog:max(t),k),2*Nf2);
     spectre_fog=spectre_fog(:);
     spec(:,k)=spec(:,k)+abs(spectre_fog).^2; 
    end;    
    % sp=abs(fft(sig(min(t):max(t),k))).^2; 
    % fr1=(0.5*(0:Lt-1)/Lt)*fs;
    % fr2=(0.5*(0:2*Nf2-1)/2/Nf2)*fs;
    % spec(:,k)=interp1(fr1,sp,fr2);
   end
   spec1=sum(spec(1:Nf2,k));
   spec2=sum(spec(Nf2+1:2*Nf2,k));
   if spec2>0.1*spec1,
    isana=0;
    if ~isreal(sig(min(t):max(t),k)),
     alpha=1;
    end
   end
  end
 end

 if layout==1,                          % Time-Frequency Representation only
  set(axtfr,'Position',[0.10 0.10 0.80 0.80]);

 elseif layout==2,			% TFR + Signal
  set(axtfr,'Position',[0.10 0.10 0.80 0.55]);
  set(axsig,'Position',[0.10 0.73 0.80 0.20]); 
  axes(axsig);  
  if sigenveloppe,
   plot((min(t):max(t))/fs,real(sig(min(t):max(t),:)),...
        (min(t):max(t))/fs, abs(sig(min(t):max(t),:)));
  else
   plot((min(t):max(t))/fs,real(sig(min(t):max(t),:)));
  end;
  set(gca,['X' TickLabelStr],[]);
  ylabel('Real part');
  title('Signal in time');
  Min=min(min(real(sig))); Max=max(max(real(sig)));
  axis([min(ts) max(ts) Min Max]);

 elseif layout==3,			% TFR + spectrum
  set(axspec,'Position',[0.10 0.10 0.15 0.80]);
  set(axtfr ,'Position',[0.35 0.10 0.55 0.80]);

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
  if linlogspec==0,
   plot(freqs,spec);
   title('Linear scale'); 
   % set(axspec,'ytick',[0 Maxsp*max(eps,threshold)*0.01 Maxsp]);
   set(axspec,'YTickMode', 'auto');
   set(axspec,'Ylim', [Maxsp*threshold*0.01 Maxsp*1.2]);
   set(axspec,'Xlim', [fmin fmax]);
  else
   plot(freqs,10*log10(spec/Maxsp)); 
   title('Log. scale [dB]');
   set(axspec,'YTickMode', 'auto');
   set(axspec,'Ylim',[10*log10(threshold*0.01) 0]);
   set(axspec,'Xlim', [fmin fmax]);
  end
  xlabel('Energy spectral density');
  Nsp=length(spec); 
  set(gca, ['X' TickLabelStr],[],'view',[-90 90]);

 elseif layout==4,			% TFR + signal + spectrum

  set(axspec,'Position',[0.10 0.10 0.15 0.55]);
  set(axsig ,'Position',[0.35 0.73 0.55 0.20]);
  set(axtfr ,'Position',[0.35 0.10 0.55 0.55]);

  axes(axsig); 
  if sigenveloppe,
   plot((min(t):max(t))/fs,real(sig(min(t):max(t),:)),...
        (min(t):max(t))/fs, abs(sig(min(t):max(t),:)));
  else
   plot((min(t):max(t))/fs,real(sig(min(t):max(t),:)));
  end;
  ylabel('Real part');
  title('Signal in time');
  set(gca,['X' TickLabelStr],[]);

  Min=min(min(real(sig))); Max=max(max(real(sig)));
  axis([min(ts) max(ts) Min Max]);

  axes(axspec); 
  if isaffine, 
% IF YOU UNDERSTAND SOMETHING TO THAT, PLEASE EXPLAIN ME (f. auger)
%   f1=freqr(1); f2=freqr(Nf2); df=f2-f1; 
%   Nf4=round((Nf2-1)*fs/(2*df))+1;
%   for k=1:Ncol,
%    spec(1:alpha*Nf4,k)=abs(fft(sig(min(t):max(t),k),alpha*Nf4)).^2; 
%   end
%   spec=spec((round(f1*2*(Nf4-1)/fs)+1):(round(f1*2*(Nf4-1)/fs)+Nf2),:);
%   freqs=linspace(f1,f2,Nf2);
   for k=1:Ncol,
    freqs=linspace(freqr(1),freqr(Nf2),Nf2); 
    spec=interp1(0.5*fs*(0:Nf2-1)/Nf2,spec(1:Nf2,k),freqs);
   end;
  else
   freqs=freqr;
   spec=spec(1:Nf2,:);
  end
  Maxsp=max(max(spec));
  if linlogspec==0,
   plot(freqs,spec);
   title('Linear scale');  
   set(axspec,'YTickMode', 'auto');
   set(axspec,'Ylim', [Maxsp*threshold*0.01 Maxsp*1.2]);
   set(axspec,'Xlim', [fmin*fs fmax*fs]);
  else
   plot(freqs,10*log10(spec/Maxsp)); 
   title('Log. scale [dB]');
   set(axspec,'Ytickmode','auto');
   set(axspec,'Ylim',[10*log10(threshold*0.01) 0]);
   set(axspec,'Xlim', [fmin*fs fmax*fs]);
  end
  xlabel('Energy spectral density');
  Nsp=length(spec); 
  set(gca,['X' TickLabelStr],[],'view',[-90 90]);

 end;

  if iscolorbar,                     % Is there a color bar ?
   PositionTfr=get(axtfr,'Position');
   set(axtfr,'Position',PositionTfr-[0 0 0.03 0]);
   set(axcb, 'Position',[PositionTfr(1)+PositionTfr(3)-0.01,...
                         PositionTfr(2) 0.01 PositionTfr(4)]);
   axes(axcb); 
   Ncolors=length(colormap); 
   [cmin,cmax]=caxis; colorvect=linspace(cmax,cmin,Ncolors);
   imagesc(colorvect'); axis('off'); 
   if issig,                          % there is a signal
    PositionSig=get(axsig,'Position');
    set(axsig,'Position',PositionSig-[0 0 0.03 0]);
   end 
  end

 axes(axtfr);                         % Display the tfr
 if (tcol==1),
  plot(freqr,tfr(1:Nf2));
  set(axtfr,'Xlim', [fmin*fs fmax*fs]);

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
  if isaffine & (display==2),
   fprintf('imagesc does not support non-linear scales for axes. Replaced by pcolor.\n');
   display=3; 
  end;

  if display==1,                  % contour
   contour(ts,freqr,tfr,levels);  % contour(tfr,levels,ts,freqr);
   set(axtfr,'Ylim', [fmin*fs fmax*fs]);
   DisplayStr=', contour';
  elseif display==2,              % imagesc
   if linlogtfr==0,
    imagesc(ts,freqr,tfr); axis('xy');
   else
    imagesc(ts,freqr,log10(tfr));axis('xy');
   end
   set(axtfr,'Ylim', [fmin*fs fmax*fs]);
   DisplayStr=', imagesc';

  elseif display==3,              % pcolor
   if linlogtfr==0,
    pcolor(ts,freqr,tfr); shading interp;
   else
    pcolor(ts,freqr,log10(tfr)); shading interp;
   end
   set(axtfr,'Ylim', [fmin*fs fmax*fs]);
   DisplayStr=', pcolor';

  elseif display==4,              % surf
   if linlogtfr==0,
    surf(ts,freqr,tfr); shading interp;
    zlabel('Amplitude');
    axis([ts(1) ts(tcol) fmin*fs fmax*fs mini maxi]);
   else
    surf(ts,freqr,log10(tfr)); shading interp;
    zlabel('Positive values');
    axis([ts(1) ts(tcol) fmin fmax log10(mini) log10(maxi)]);
   end
   DisplayStr=', surf';

  elseif display==5,              % mesh
   if linlogtfr==0,
    mesh(ts,freqr,tfr); shading interp;
    zlabel('Amplitude');
    axis([ts(1) ts(tcol) fmin*fs fmax*fs mini maxi]);
   else
    mesh(ts,freqr,log10(tfr));shading interp;
    zlabel('Positive values');
    axis([ts(1) ts(tcol) fmin*fs fmax*fs log10(mini) log10(maxi)]);
   end
   DisplayStr=', mesh';
  end
 
 % Define the title and check the input arguments depending on 'method'

 method=method(4:length(method));

 if nargin==5, % if there is no additional parameters, do the best.
  title([method, LinLogStr,DisplayStr,...
         ', Threshold=',num2str(threshold),'%']);

 elseif strcmp(method,'WV'  ) | strcmp(method,'MH') | ...
        strcmp(method,'PAGE'), % no parameters
  title([method,', Nf=',num2str(Nf2), LinLogStr,DisplayStr,...
        ', Threshold=',num2str(threshold),'%']);

 elseif strcmp(method,'PWV'  )|strcmp(method,'PMH'  )| ...
        strcmp(method,'SP'   )|strcmp(method,'PPAGE')| ...
        strcmp(method,'RSP'  )|strcmp(method,'RPPAG')| ...
        strcmp(method,'RPWV' )|strcmp(method,'RPMH' ),
  h=p1;[hrow,hcol]=size(h); Lh=(hrow-1)/2; % one parameter
  if (hcol~=1)|(rem(hrow,2)==0),
   error('h must be a smoothing window with odd length'); end;
  title([method, ', Lh=',num2str(Lh), ', Nf=',num2str(Nf2),...
        LinLogStr, DisplayStr,', Threshold=',num2str(threshold),'%']);
 
 elseif strcmp(method,'STFT'), % short-time fourier transform case
  h=p1;[hrow,hcol]=size(h); Lh=(hrow-1)/2;
  if (hcol~=1)|(rem(hrow,2)==0),
   error('h must be a smoothing window with odd length'); end;
  title(['|',method, '|^2, Lh=',num2str(Lh),...
        ', Nf=',num2str(Nf2), LinLogStr, DisplayStr,', Thld=',...
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
         ', Nf=',num2str(Nf2),LinLogStr, DisplayStr,...
         ', Threshold=',num2str(threshold),'%']);
 
 elseif strcmp(method,'MMCE'),
  h=p1;[hrow,hcol]=size(h); Lh=(hrow-1)/2;
  if (rem(hrow,2)==0),
   error('h must be a smoothing window with odd length'); end;
  title([method, ', Lh=',num2str(Lh), ', Nf=',num2str(Nf2),...
        LinLogStr, DisplayStr,', Threshold=',num2str(threshold),'%']);
 
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
         LinLogStr, DisplayStr, ', Threshold=',num2str(threshold),'%']);
 
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
         LinLogStr, DisplayStr, ', Threshold=',num2str(threshold),'%']);
 
 elseif strcmp(method,'MSC' ) | strcmp(method,'RMSC' )
  f0T=p1; if (f0T<=0), error('f0T must be positive'); end;
  title([method, ', f0T=',num2str(f0T), ', Nf=',num2str(Nf2),...
        LinLogStr, DisplayStr, ', Threshold=',num2str(threshold),'%']);
 
 elseif strcmp(method,'RGAB' )
  Nh=p1; if (Nh<=0), error('Nh must be positive'); end;
  title([method, ', Nh=',num2str(Nh), ', Nf=',num2str(Nf2),...
        LinLogStr, DisplayStr, ', Threshold=',num2str(threshold),'%']);
 
 
 elseif strcmp(method,'DFLA' ) | strcmp(method,'UNTER' )| ...
        strcmp(method,'BERT' ),
  N=p1; 
  if (N<=0),                error('N must be positive'); end;
  title([method, ', N=',num2str(N), LinLogStr, DisplayStr, ', Threshold=',...
         num2str(threshold), '%']);  
 
 elseif strcmp(method,'SCALO'),
  Nh0=p1; N=p2; 
  if (Nh0<0),               error('Nh0 must be positive'); end;
  if (N<=0),                error('N must be positive'); end;
  if (Nh0>0), 
   title([method, ', Morlet wavelet, Nh0=', num2str(Nh0), ...
         ', N=',num2str(N), LinLogStr, DisplayStr, ', Thld=',...
         num2str(threshold), '%']);  
  else 
   title([method, ', Mexican hat, N=',num2str(N), LinLogStr, DisplayStr, ...
          ', Thld=', num2str(threshold), '%']);  
  end
 
 elseif strcmp(method,'ASPW'),
  Nh0=p1; Ng0=p2; N=p3; 
  if (Nh0<0),               error('Nh0 must be positive'); end;
  if (Ng0<0),               error('Ng0 must be positive'); end;
  if (N<=0),                error('N must be positive'); end;
  if (Nh0>0), 
   title([method, ', Morlet wlt, Nh0=', num2str(Nh0), ', Ng0=',...
         num2str(Ng0), ', N=',num2str(N), LinLogStr, DisplayStr, ', Thld=',...
         num2str(threshold), '%']);  
  else 
   title([method, ', Mexican hat, Ng0=',num2str(Ng0),...  
         ', N=',num2str(N), LinLogStr, DisplayStr, ', Thld=',... 
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
         num2str(Ng0), ', N=',num2str(N), LinLogStr, DisplayStr, ', Thld=',...
         num2str(threshold), '%']);  
  else 
   title([method, ', K=', num2str(K), ', Mexican hat, Ng0=',...
         num2str(Ng0),', N=',num2str(N), LinLogStr, DisplayStr, ', Thld=',...
         num2str(threshold), '%']);  
  end
 
 
 elseif strcmp(method,'GABOR'),
  N=p1; Q=p2; h=p3; [hrow,hcol]=size(h); Lh=(hrow-1)/2;
  if (hcol~=1)|(rem(hrow,2)==0),
   error('h must be a smoothing window with odd length'); end;
  title([method, ', Lh=',num2str(Lh), ', Nf=',...
         num2str(Nf2),', N=',num2str(N),', Q=',num2str(Q),...
         LinLogStr, DisplayStr, ', Thld=',num2str(threshold),'%']);
 
 end;
end

% add the correct legend on the axes
if unitHz==1,
 xlabel('Time [s]'); ylabel('Frequency [Hz]');
elseif unitHz==2,
 xlabel('Time [ms]'); ylabel('Frequency [kHz]');
elseif unitHz==3,
 xlabel('Time [µs]'); ylabel('Frequency [MHz]');
end


if (isgridsig & issig),		% Updating of the grids
 axes(axsig); grid on	
elseif (~isgridsig & issig),
 axes(axsig); grid off
end 
if (isgridspec & isspec),
 axes(axspec); grid on
elseif (~isgridspec & isspec),
 axes(axspec); grid off
end 

if (isgridtfr),			% upating of the grid on the tfr
 axes(axtfr); grid on  
elseif (~isgridtfr),
 axes(axtfr); grid off
end 

