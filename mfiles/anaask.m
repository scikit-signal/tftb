function [y,am]=anaask(N,Ncomp,f0);
%ANAASK	Amplitude Shift Keying (ASK) signal.
% 	[Y,AM]=ANAASK(N,NCOMP,F0) returns a complex amplitude
%	modulated signal of normalized frequency F0, with a uniformly
%	distributed random amplitude. 
% 	Such signal is only 'quasi'-analytic.
%	
%	N     : number of points
%	NCOMP : number of points of each component (default: N/5)
% 	F0    : normalized frequency.              (default: 0.25)
% 	Y     : signal
% 	AM    : resulting amplitude modulation     (optional).
%
%	Example :
%	 [signal,am]=anaask(512,64,0.05); clf; figure(gcf);
%  	 subplot(211); plot(real(signal)); subplot(212); plot(am);
%
%	See also ANAFSK, ANABPSK, ANAQPSK.

%	O. Lemoine - October 1995
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
 error('The number of parameters must be at least 1.');
elseif (nargin == 1),
 Ncomp=round(N/5); f0=0.25;
elseif (nargin == 2),
 f0=0.25;
end;

if (N <= 0),
 error('The signal length N must be strictly positive' );
elseif (f0<0)|(f0>0.5),
 error('f0 must be between 0 and 0.5');
end;

MatlabVersion=version; MatlabVersion=str2num(MatlabVersion(1));
if (MatlabVersion==4), rand('uniform'); end;

m=ceil(N/Ncomp); jumps=rand(m,1);
am=kron(jumps,ones(Ncomp,1)); am=am(1:N,1);
y=am.*fmconst(N,f0,1);
