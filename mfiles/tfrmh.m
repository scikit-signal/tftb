function [tfr,t,f] = tfrmh(x,t,N,trace);
%TFRMH	Margenau-Hill time-frequency distribution.
%	[TFR,T,F]=TFRMH(X,T,N,TRACE) computes the Margenau-Hill
%	distribution of a discrete-time signal X, 
%	or the cross Margenau-Hill representation between two signals. 
% 
%	X     : signal if auto-MH, or [X1,X2] if cross-MH.
%	T     : time instant(s)          (default : 1:length(X)).
%	N     : number of frequency bins (default : length(X)).
%	TRACE : if nonzero, the progression of the algorithm is shown
%                                        (default : 0).
%	TFR   : time-frequency representation. When called without 
%               output arguments, TFRMH runs TFRQVIEW.
%	F     : vector of normalized frequencies.
%
%	Example :
%	 sig=fmlin(128,0.1,0.4); tfrmh(sig,1:128,128,1);
% 
%	See also  all the time-frequency representations listed in
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
if (nargin == 1),
 t=1:xrow; N=xrow; trace=0;
elseif (nargin == 2),
 N=xrow; trace=0;
elseif (nargin == 3),
 trace = 0;
end;

if (N<0),
 error('N must be greater than zero');
end;
[trow,tcol] = size(t);
if (xcol==0) | (xcol>2),
 error('X must have one or two columns');
elseif (trow~=1),
 error('T must only have one row'); 
elseif (2^nextpow2(N)~=N),
 fprintf('For a faster computation, N should be a power of two\n');
end; 

tfr= zeros (N,tcol) ;  
if trace, disp('Margenau-Hill distribution'); end;
for icol=1:tcol,
 ti= t(icol); tau=-min([N-ti,xrow-ti]):(ti-1);
 indices= rem(N+tau,N)+1; 
 if trace, disprog(icol,tcol,10); end;
 tfr(indices,icol)=x(ti,1)*conj(x(ti-tau,xcol));
end; 

tfr= real(fft(tfr)); 

if (nargout==0),
 tfrqview(tfr,x,t,'tfrmh');
elseif (nargout==3),
 if rem(N,2)==0, 
  f=[0:N/2-1 -N/2:-1]'/N;
 else
  f=[0:(N-1)/2 -(N-1)/2:-1]'/N;  
 end;
end;

