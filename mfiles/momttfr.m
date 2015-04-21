function [fm,B2]=momttfr(tfr,method,fbmin,fbmax,freqs);
%MOMTTFR Time moments of a time-frequency representation.
%	[FM,B2]=MOMTTFR(TFR,METHOD,FBMIN,FBMAX,FREQS) computes the time
%	moments of a time-frequency representation.
% 
%	TFR   : time-frequency representation ([Nrow,Ncol]=size(TFR)).
%	METHOD: chosen representation (name of the corresponding M-file).  
%	FBMIN : smallest frequency bin (default : 1)
%	FBMAX : highest  frequency bin (default : Nrow)
%	FREQS : true frequency of each frequency bin. FREQS must be of
%		length FBMAX-FBMIN+1. 
%	         (default : 0:step:(0.5-step) or -0.5:step:(0.5-step) 
%		             depending on METHOD) 
%	FM    : averaged frequency     (first order moment)
%	B2    : frequency band         (second order moment)
%
%	Example :
%	 sig=fmlin(128,0.1,0.4); tfr=tfrwv(sig);
%	 [fm,B2]=momttfr(tfr,'tfrwv'); 
%	 subplot(211); plot(fm); subplot(212); plot(B2);
%	 freqs=linspace(0,63/128,64); tfr=tfrsp(sig); 
%	 [fm,B2]=momttfr(tfr,'tfrsp',1,64,freqs); 
%	 subplot(211); plot(fm); subplot(212); plot(B2);
%
%	See also MOMFTFR, MARGTFR.

%	F. Auger, August 1995.
%	Copyright (c) 1996 by CNRS (France).
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

[tfrrow,tfrcol]=size(tfr);
if (nargin<=1),
 error('At least two input arguments');
elseif (nargin==2),
 fbmin=1; fbmax=tfrrow; 
elseif (nargin==3),
 fbmax=tfrrow; 
end
if (fbmin>fbmax)|(fbmin==0)|(fbmax>tfrrow),
 error('1<=FBMIN<=FBMAX<=Nrow');
elseif nargin==5,
 if length(freqs)~=(fbmax-fbmin+1),
  error('FREQS must have FBMAX-FBMIN+1 elements');
 end
end;

method=upper(method);
if strcmp(method,'TFRWV'   ) | strcmp(method,'TFRPWV'  ) | ...
   strcmp(method,'TFRSPWV' ) | strcmp(method,'TFRCW'   ) | ...
   strcmp(method,'TFRZAM'  ) | strcmp(method,'TFRBJ'   ) | ...
   strcmp(method,'TFRBUD'  ) | strcmp(method,'TFRGRD'  ) | ...
   strcmp(method,'TFRRSPWV') | strcmp(method,'TFRRPWV' ) | ...
   strcmp(method,'TFRRIDB' ) | strcmp(method,'TFRRIDH' ) | ...
   strcmp(method,'TFRRIDT' ) | strcmp(method,'TFRASPW' ) | ...
   strcmp(method,'TFRDFLA' ) | strcmp(method,'TFRSPBK' ) | ...
   strcmp(method,'TFRUNTAC') | strcmp(method,'TFRUNTPA') | ...
   strcmp(method,'TFRBERT' ) | strcmp(method,'TFRSCALO') | ...
   strcmp(method,'TYPE2' ),
 typertf='TYPE2';
elseif strcmp(method,'TFRPMH'  )| strcmp(method,'TFRRPMH' )| ...
       strcmp(method,'TFRSP'   )| strcmp(method,'TFRRSP'  )| ...
       strcmp(method,'TFRPPAGE')| strcmp(method,'TFRRPPAG')| ...
       strcmp(method,'TFRMHS'  )| strcmp(method,'TFRRGAB' )| ...
       strcmp(method,'TFRMH'   )| strcmp(method,'TFRMMCE' )| ...
       strcmp(method,'TFRMSC'  )| strcmp(method,'TFRRMSC' )| ...
       strcmp(method,'TFRPAGE' )| strcmp(method,'TFRGABOR')| ...
       strcmp(method,'TFRRI'   )| strcmp(method,'TYPE1'   ),
 typertf='TYPE1';
else
 error('Unknown representation.');
end;

if nargin<=4, 
 if strcmp(typertf,'TYPE1'),
  freqs=rem((fbmin-1:fbmax-1)/tfrrow + 0.5, 1.0) - 0.5;
 elseif strcmp(typertf,'TYPE2'),
  freqs=0.5*(fbmin-1:fbmax-1)/tfrrow;
 end;
end

E  = sum(tfr(fbmin:fbmax,:));
fm = (freqs * tfr(fbmin:fbmax,:) ./E).'; 
B2 = (freqs.^2 * tfr(fbmin:fbmax,:) ./E).' - fm.^2;


