function [tfr,t,f,wt]=tfrscalo(X,time,wave,fmin,fmax,N,trace);
%TFRSCALO Scalogram, for Morlet or Mexican hat wavelet.
%	[TFR,T,F,WT]=TFRSCALO(X,T,WAVE,FMIN,FMAX,N,TRACE) computes 
%	the scalogram (squared magnitude of a continuous wavelet
%	transform). 
%
%	X : signal (in time) to be analyzed (Nx=length(X)). Its
%	    analytic version is used (z=hilbert(real(X))).  
%	T : time instant(s) on which the TFR is evaluated 
%	     					(default : 1:Nx).
%	WAVE : half length of the Morlet analyzing wavelet at coarsest 
% 	    scale. If WAVE = 0, the Mexican hat is used. WAVE can also be
%           a vector containing the time samples of any bandpass
%           function, at any scale.        	(default : sqrt(Nx)). 
%	FMIN,FMAX : respectively lower and upper frequency bounds of 
%	    the analyzed signal. These parameters fix the equivalent
%	    frequency bandwidth (expressed in Hz). When unspecified, you
%	    have to enter them at the command line from the plot of the
%	    spectrum. FMIN and FMAX must be >0 and <=0.5.
%	N : number of analyzed voices.
%	TRACE : if nonzero, the progression of the algorithm is shown
%	                                 	(default : 0).
%	TFR : time-frequency matrix containing the coefficients of the
%	    decomposition (abscissa correspond to uniformly sampled time,
%	    and ordinates correspond to a geometrically sampled
%	    frequency). First row of TFR corresponds to the lowest 
%	    frequency. When called without output arguments, TFRSCALO
%	    runs TFRQVIEW.
%	F : vector of normalized frequencies (geometrically sampled 
%	    from FMIN to FMAX).
%	WT : Complex matrix containing the corresponding wavelet
%	    transform. The scalogram TFR is the square modulus of WT.
%
%	Example :    
%	 sig=altes(64,0.1,0.45); tfrscalo(sig);  
%
%	See also all the time-frequency representations listed in
%	the file CONTENTS (TFR*)

%	P. Goncalves, October 1995 - O. Lemoine, June 1996. 
%	Copyright (c) 1995 Rice University - CNRS 1996.
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

if (nargin == 0),
 error('At least one parameter required');
end;

[xrow,xcol] = size(X);
if nargin<=6, trace=0; end

if (nargin == 1),
 time=1:xrow; wave=sqrt(xrow);
elseif (nargin == 2),
 wave=sqrt(xrow);
elseif (nargin==4),
 disp('FMIN will not be taken into account. Determine it with FMAX');
 disp('     from the following plot of the spectrum.'); 
elseif nargin==5,
 N=xrow;
end;

[trow,tcol] = size(time);
if (xcol==0)|(xcol>2),
 error('X must have one or two columns');
elseif (trow~=1),
 error('TIME must only have one row'); 
elseif wave<0,
 error('WAVE must be positive');
end; 

s = (real(X) - mean(real(X)))';  
z = hilbert(s) ;

if trace, disp('Scalogram distribution'); end;

if nargin<=4		        % fmin,fmax,N unspecified
 STF = fft(fftshift(z(min(time):max(time)))); Nstf=length(STF);
 sp = (abs(STF(1:round(Nstf/2)))).^2; Maxsp=max(sp);
 f = linspace(0,0.5,round(Nstf/2)+1) ; f = f(1:round(Nstf/2));
 plot(f,sp) ; grid;
 xlabel('Normalized frequency');
 title('Analyzed signal energy spectrum');
 indmin=min(find(sp>Maxsp/100));
 indmax=max(find(sp>Maxsp/100));
 fmindflt=max([0.01 0.05*fix(f(indmin)/0.05)]);
 fmaxdflt=0.05*ceil(f(indmax)/0.05);
 txtmin=['Lower frequency bound [',num2str(fmindflt),'] : '];
 txtmax=['Upper frequency bound [',num2str(fmaxdflt),'] : '];
 fmin = input(txtmin); fmax = input(txtmax);
 if isempty(fmin), fmin=fmindflt; end
 if isempty(fmax), fmax=fmaxdflt; end
 txt=['Number of frequency samples [',num2str(2^nextpow2(xrow)),'] : ']; 
 N=input(txt); 
 if isempty(N), N=2^nextpow2(xrow); end
end

fmin_s=num2str(fmin); fmax_s=num2str(fmax); 
N_s=num2str(N);

if fmin >= fmax
 error('FMAX must be greater or equal to FMIN');
elseif fmin<=0.0 | fmin>0.5,
 error('FMIN must be > 0 and <= 0.5');
elseif fmax<=0.0 | fmax>0.5,
 error('FMAX must be > 0 and <= 0.5');
end
if trace,
 disp(['Frequency runs from ',fmin_s,' to ',fmax_s,' with ',N_s,' points']);
end

f = logspace(log10(fmin),log10(fmax),N);
a = logspace(log10(fmax/fmin),log10(1),N); 


wt =zeros(N,tcol);
tfr=zeros(N,tcol);

if wave > 0
 if trace, disp(['using a Morlet wavelet']) ; end
 for ptr=1:N,
  if trace, disprog(ptr,N,10); end
  nha = wave*a(ptr);
  tha = -round(nha) : round(nha);
  ha  = exp(-(2*log(10)/nha^2)*tha.^2).*exp(i*2*pi*f(ptr)*tha); 
  detail = conv(z,ha)./sqrt(a(ptr));
  detail = detail(round(nha)+1:length(detail)-round(nha)) ;
  wt(ptr,:)  = detail(time) ;
  tfr(ptr,:) = detail(time).*conj(detail(time)) ;
 end
elseif wave == 0
 if trace, disp(['using a Mexican hat wavelet']) ; end
 for ptr = 1:N
  if trace, disprog(ptr,N,10); end
  ha  = mexhat(f(ptr)) ;
  nha = (length(ha)-1)/2 ;
  detail = conv(z,ha)./sqrt(a(ptr));
  detail = detail(round(nha)+1:length(detail)-round(nha)) ;
  wt(ptr,:)  = detail(time);
  tfr(ptr,:) = detail(time).*conj(detail(time)) ;
 end  
elseif length(wave) > 1
 [rwav,cwav]=size(wave);
 if cwav>rwav, wave=wave.'; end
 wavef = fft(wave) ;
 nwave = length(wave) ;
 f0 = find(abs(wavef(1:nwave/2)) == max(abs(wavef(1:nwave/2))));
 f0 = mean((f0-1).*(1/nwave));
 if trace, disp(['mother wavelet centered at f0 = ',num2str(f0)]); end
 a = logspace(log10(f0/fmin),log10(f0/fmax),N);
 B = 0.99;
 R = B/((1.001)/2); 
 nscale = max(128,round((B*nwave*(1+2/R)*log((1+R/2)/(1-R/2)))/2));
 if trace, disp('Scale computation :'); end
 wts = scale(wave,a,fmin,fmax,nscale,trace);
 for ptr = 1:N, 
  clear detail
  if trace, disprog(ptr,N,10); end
  ha = wts(ptr,:);
  nha = length(ha)/2;
  detail = conv(z,ha)./sqrt(a(ptr));
  detail = detail(fix(nha):length(detail)-round(nha));
  wt(ptr,:) = detail(time);
  tfr(ptr,:) = detail(time).*conj(detail(time));
 end
end


t = time;
f = f';

% Normalization
SP = fft(z); 
indmin = 1+round(fmin*(xrow-2));
indmax = 1+round(fmax*(xrow-2));
SPana=SP(indmin:indmax);
tfr=tfr*norm(SPana)^2/integ2d(tfr,t,f)/N;

if (nargout==0),
 tfrqview(tfr,hilbert(real(X)),t,'tfrscalo',wave,N,f);
end;


