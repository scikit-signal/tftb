function [tfr,t,f] = tfrmhs(x,t,N,g,h,trace);
%TFRMHS	Margenau-Hill-Spectrogram time-frequency distribution.
%	[TFR,T,F]=TFRMHS(X,T,N,G,H,TRACE) computes the Margenau-Hill-Spectrogram 
%	distribution of a discrete-time signal X, or the cross
%	Margenau-Hill-Spectrogram representation between two signals. 
% 
%	X     : Signal if auto-MHS, or [X1,X2] if cross-MHS.
%	T     : time instant(s)          (default : 1:length(X)).
%	N     : number of frequency bins (default : length(X)).
%	G,H   : analysis windows, normalized so that the representation 
%               preserves the signal energy.
%	                (default : Hamming(N/10) and Hamming(N/4)). 
%	TRACE : if nonzero, the progression of the algorithm is shown
%                                        (default : 0).
%	TFR   : time-frequency representation. When called without 
%               output arguments, TFRMHS runs TFRQVIEW.
%	F     : vector of normalized frequencies.
%
%	Example:
%	 sig=fmlin(128,0.1,0.4); g=tftb_window(21,'Kaiser'); 
%	 h=tftb_window(63,'Kaiser'); tfrmhs(sig,1:128,64,g,h,1);
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

[grow,gcol]=size(g); Lg=(grow-1)/2;
if (gcol~=1)|(rem(grow,2)==0),
 error('G must be a smoothing window with odd length'); 
end;

[hrow,hcol]=size(h); Lh=(hrow-1)/2; h=h/h(Lh+1);
if (hcol~=1)|(rem(hrow,2)==0),
  error('H must be a smoothing window with odd length');
end;

Lgh=min(Lg,Lh); points=-Lgh:Lgh; 
Kgh=sum(h(Lh+1+points).*conj(g(Lg+1+points))); h=h/Kgh;

tfr= zeros (N,tcol); tfr2= zeros(N,tcol);
if trace, disp('Pseudo Margenau-Hill distribution'); end;
for icol=1:tcol,
 ti= t(icol);
 if trace, disprog(icol,tcol,10); end;
 tau=-min([round(N/2)-1,Lg,ti-1]):min([round(N/2)-1,Lg,xrow-ti]);
 indices= rem(N+tau,N)+1;
 tfr(indices,icol)=x(ti+tau,1).*conj(g(Lg+1+tau));
 tau=-min([round(N/2)-1,Lh,ti-1]):min([round(N/2)-1,Lh,xrow-ti]);
 indices= rem(N+tau,N)+1;
 tfr2(indices,icol)=x(ti+tau,xcol).*conj(h(Lh+1+tau));
end; 
if trace, fprintf('\n'); end;
tfr=real(fft(tfr).*conj(fft(tfr2))); 

if (nargout==0),
 tfrqview(tfr,x,t,'tfrmhs',g,h);
elseif (nargout==3),
 if rem(N,2)==0, 
  f=[0:N/2-1 -N/2:-1]'/N;
 else
  f=[0:(N-1)/2 -(N-1)/2:-1]'/N;  
 end;
end;


