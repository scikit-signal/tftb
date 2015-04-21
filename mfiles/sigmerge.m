function sig=sigmerge(x1,x2,ratio);
%SIGMERGE Add two signals with given energy ratio in dB.
%	SIG=SIGMERGE(X1,X2,RATIO) adds two signals so that a given
%	energy ratio expressed in deciBels is satisfied.
% 
%	X1, X2 : input signals.
%	RATIO  : Energy ratio in deciBels	(default : 0 dB).
%	X      : output signal.
%	X= X1+H*X2, such that 10*log(Energy(X1)/Energy(H*X2))=RATIO
%
%	Example : 
%	 sig=fmlin(64,0.01,0.05,1); noise=hilbert(randn(64,1));
%	 SNR=15; x=sigmerge(sig,noise,SNR);
%	 Esig=mean(abs(sig).^2); Enoise=mean(abs(x-sig).^2);
%	 10*log10(Esig/Enoise)

%	F. Auger, July 1995.
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

if (nargin<2)
 error('At least two parameters are required');
elseif nargin==2,
 ratio=0;
end;

[x1row,x1col] = size(x1);
[x2row,x2col] = size(x2);

if (x1col~=1)|(x2col~=1),
 error('X1 and X2 must have only one column');
elseif (x1row~=x2row),
 error('X1 and X2 must have the same number of rows');
elseif (length(ratio)~=1),
 error('RATIO must be a scalar');
elseif (ratio==inf),
 sig = x1;
else
 Ex1=mean(abs(x1).^2);
 Ex2=mean(abs(x2).^2);
 h=sqrt(Ex1/(Ex2*10^(ratio/10)));
 sig=x1+h*x2;
end;
