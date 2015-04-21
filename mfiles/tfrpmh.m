function [tfr,t,f] = tfrpmh(x,t,N,h,trace);
%TFRPMH  Pseudo Margenau-Hill time-frequency distribution.
%       [TFR,T,F]=TFRPMH(X,T,N,H,TRACE) computes the Pseudo Margenau-Hill  
%       distribution of a discrete-time signal X, or the
%       cross Pseudo Margenau-Hill representation between two signals. 
% 
%       X     : signal if auto-PMH, or [X1,X2] if cross-PMH.
%       T     : time instant(s)          (default : 1:length(X)).
%       N     : number of frequency bins (default : length(X)).
%       H     : frequency smoothing window, H(0) being forced to 1
%                                        (default : Hamming(N/4)). 
%       TRACE : if nonzero, the progression of the algorithm is shown
%                                        (default : 0).
%       TFR   : time-frequency representation. When called without 
%               output arguments, TFRPMH runs TFRQVIEW.
%       F     : vector of normalized frequencies.
%
%       Example :
%        sig=fmlin(128,0.1,0.4); t=1:128; 
%        h=tftb_window(63,'Kaiser'); tfrpmh(sig,t,128,h,1);
% 
%       See also all the time-frequency representations listed in
%        the file CONTENTS (TFR*)

%       F. Auger, May-August 1994, July 1995.
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

[xrow,xcol] = size(x);
if (nargin < 1),
 error('At least 1 parameter is required');
elseif nargin<=2,
 N=xrow;
end;

hlength=floor(N/4);
hlength=hlength+1-rem(hlength,2);

if (nargin == 1),
 t=1:xrow; h = tftb_window(hlength); trace=0;
elseif (nargin == 2 | nargin == 3),
 h = tftb_window(hlength); trace = 0;
elseif (nargin == 4),
 trace = 0;
end;

if (N<0),
 error('N must be greater than zero');
end;
[trow,tcol] = size(t);
if (xcol==0)|(xcol>2),
 error('X must have one or two columns');
elseif (trow~=1),
 error('T must only have one row'); 
elseif (2^nextpow2(N)~=N), %(rem(log(N)/log(2),1)~=0),
 fprintf('For a faster computation, N should be a power of two\n');
end; 

[hrow,hcol]=size(h); Lh=(hrow-1)/2; h=h/h(Lh+1);
if (hcol~=1)|(rem(hrow,2)==0),
 error('H must be a smoothing window with odd length');
end;

tfr= zeros (N,tcol) ;  
if trace, disp('Pseudo Margenau-Hill distribution'); end;
for icol=1:tcol,
 ti= t(icol); tau=-min([round(N/2)-1,Lh,xrow-ti]):min([round(N/2)-1,Lh,ti-1]);
 indices= rem(N+tau,N)+1;
 if trace, disprog(icol,tcol,10); end;
 tfr(indices,icol)=h(Lh+1+tau).*x(ti,1).*conj(x(ti-tau,xcol));
end; 
if trace, fprintf('\n'); end;
tfr= real(fft(tfr)); 

if (nargout==0),
 tfrqview(tfr,x,t,'tfrpmh',h);
elseif (nargout==3),
 if rem(N,2)==0, 
  f=[0:N/2-1 -N/2:-1]'/N;
 else
  f=[0:(N-1)/2 -(N-1)/2:-1]'/N;  
 end;
end;

