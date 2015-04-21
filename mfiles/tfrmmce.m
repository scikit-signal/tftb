function [tfr,t,f] = tfrmmce(x,h,t,N,trace);
%TFRMMCE Minimum mean cross-entropy combination of spectrograms.
%	[TFR,T,F]=TFRMMCE(X,H,T,N,TRACE) computes the minimum mean 
%	cross-entropy combination of spectrograms using as 
%	windows the columns of the matrix H.
% 
%	X     : signal.
%	H     : frequency smoothing windows, the H(:,i) being normalized
%	        so as to be of unit energy. 
%	T     : time instant(s)          (default : 1:length(X)).
%	N     : number of frequency bins (default : length(X))
%	TRACE : if nonzero, the progression of the algorithm is shown
%                                        (default : 0).
%	TFR   : time-frequency representation. When called without 
%               output arguments, TFRMMCE runs TFRQVIEW.
%	F     : vector of normalized frequencies.
%
%	Example :
%	 sig=fmlin(128,0.1,0.4); h=zeros(19,3);
%	 h(10+(-5:5),1)=tftb_window(11); h(10+(-7:7),2)=tftb_window(15);  
%	 h(10+(-9:9),3)=tftb_window(19); tfrmmce(sig,h);
% 
%	See also all the time-frequency representations listed in
%	 the file CONTENTS (TFR*)

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

[xrow,xcol] = size(x);
if (nargin < 2),
 error('At least 2 parameters are required');
elseif (nargin==2),
 t=1:xrow; N=xrow; trace=0;
elseif (nargin==3),
 N=xrow; trace=0;
elseif (nargin==4),
 trace=0;
end;
[trow,tcol] = size(t);

if (xcol~=1),
 error('X must have only one column');
elseif (trow~=1),
 error('T must only have one row');
elseif (N<=0),
 error('N must be greater than zero'); 
elseif (2^nextpow2(N)~=N & nargin==5),
 fprintf('For a faster computation, N should be a power of two\n');
end; 

[hrow,hcol]=size(h); Lh=(hrow-1)/2; 
if (rem(hrow,2)==0),
 error('H must have an odd number of lines');
elseif hcol==1, 
 error('H must have at least two columns');
end;
h=h*diag(1.0 ./ sqrt(sum(abs(h).^2)));

tfr= zeros (N,tcol);
slides= zeros(N,hcol);
if trace, disp('MMCE Spectrogram'); end;

for icol=1:tcol,
 ti= t(icol); tau=-min([round(N/2)-1,Lh,ti-1]):min([round(N/2)-1,Lh,xrow-ti]);
 indices= rem(N+tau,N)+1; slides= zeros(N,hcol);
 if trace, disprog(icol,tcol,10); end;
 for ih=1:hcol,
  slides(indices,ih)=x(ti+tau).*conj(h(Lh+1+tau,ih));
 end;
 slides=abs(fft(slides)).^2;
 tfr(:,icol)=prod(slides')' .^(1/hcol);
end;

tfr=tfr*(sum(abs(x(t)).^2)/sum(sum(tfr)));
if trace, fprintf('\n'); end;

if (nargout==0),
 tfrqview(tfr,x,t,'tfrmmce',h);
elseif (nargout==3),
 f=fftshift((-(N/2):(N/2)-1)/N)';
end;

