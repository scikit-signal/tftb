function [tfr,t,f] = tfrzam(x,t,N,g,h,trace);
%TFRZAM	Zao-Atlas-Marks time-frequency distribution.
%	[TFR,T,F]=TFRZAM(X,T,N,G,H,TRACE) computes the Zao-Atlas-Marks 
%	distribution of a discrete-time signal X, or the
%	cross Zao-Atlas-Marks representation between two signals. 
% 
%	X     : signal if auto-ZAM, or [X1,X2] if cross-ZAM.
%	T     : time instant(s)          (default : 1:length(X)).
%	N     : number of frequency bins (default : length(X)).
%	G     : time smoothing window, G(0) being forced to 1. 
%	                                 (default : Hamming(N/10)). 
%	H     : frequency smoothing window, H(0) being forced to 1.
%	                                 (default : Hamming(N/4)). 
%	TRACE : if nonzero, the progression of the algorithm is shown
%	                                 (default : 0).
%	TFR   : time-frequency representation. When called without 
%	        outpout arguments, TFRZAM runs TFRQVIEW.
%	F     : vector of normalized frequencies.
%
%	Example :
%	 sig=fmlin(128,0.05,0.3)+fmlin(128,0.15,0.4);  
%	 g=tftb_window(9,'Kaiser'); h=tftb_window(27,'Kaiser'); 
%	 tfrzam(sig,1:128,128,g,h,1);
% 
%	See also all the time-frequency representations listed in
%	 the file CONTENTS (TFR*)

%	F. Auger, May-August 1994, July 1995.
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

if (nargin == 0),
 error('At least 1 parameter required');
end;
[xrow,xcol] = size(x);
if (xcol==0)|(xcol>2),
 error('X must have one or two columns');
end

if (nargin <= 2),
 N=xrow;
elseif (N<0),
 error('N must be greater than zero');
elseif (2^nextpow2(N)~=N),
 fprintf('For a faster computation, N should be a power of two\n');
end;

hlength=floor(N/4); hlength=hlength+1-rem(hlength,2); 
glength=floor(N/10);glength=glength+1-rem(glength,2);

if (nargin == 1),
 t=1:xrow; g = tftb_window(glength); h = tftb_window(hlength); trace = 0;
elseif (nargin == 2)|(nargin == 3),
 g = tftb_window(glength); h = tftb_window(hlength); trace = 0;
elseif (nargin == 4),
 h = tftb_window(hlength); trace = 0;
elseif (nargin == 5),
 trace = 0;
end;

[trow,tcol] = size(t);
if (trow~=1),
 error('T must only have one row'); 
end; 

[grow,gcol]=size(g); Lg=(grow-1)/2; g=g/g(Lg+1);
if (gcol~=1)|(rem(grow,2)==0),
  error('G must be a smoothing window with odd length'); 
end;

[hrow,hcol]=size(h); Lh=(hrow-1)/2; h=h/h(Lh+1);
if (hcol~=1)|(rem(hrow,2)==0),
  error('H must be a smoothing window with odd length');
end;

tfr= zeros (N,tcol) ;  
if trace, disp('Zao-Atlas-Marks Representation'); end;
for icol=1:tcol,
 ti= t(icol); taumax=min([ti+Lg-1,xrow-ti+Lg,round(N/2)-1,Lh]);
 if trace, disprog(icol,tcol,10); end;
 tfr(1,icol)= g(Lg+1) * x(ti,1) .* conj(x(ti,xcol));

 for tau=1:taumax,
  points= -min([tau,Lg,xrow-ti-tau]):min([tau,Lg,ti-tau-1]);
  g2=g(Lg+1+points);
  R=sum(g2 .* x(ti+tau-points,1) .* conj(x(ti-tau-points,xcol)));
  tfr(  1+tau,icol)=h(Lh+tau+1)*R;
  R=sum(g2 .* x(ti-tau-points,1) .* conj(x(ti+tau-points,xcol)));
  tfr(N+1-tau,icol)=h(Lh-tau+1)*R;
 end;
end; 

if trace, fprintf('\n'); end;
tfr= fft(tfr); 
if (xcol==1), tfr=real(tfr); end ;

if (nargout==0),
 tfrqview(tfr,x,t,'tfrzam',g,h);
elseif (nargout==3),
 f=(0.5*(0:N-1)/N)';
end;

