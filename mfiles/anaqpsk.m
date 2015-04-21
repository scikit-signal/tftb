function [y,pm0]=anaqpsk(N,Ncomp,f0);
%ANAQPSK Quaternary Phase Shift Keying (QPSK) signal.
% 	[Y,PM]=ANAQPSK(N,NCOMP,F0) returns a complex phase modulated signal
% 	of normalized frequency F0, whose phase changes every NCOMP point according
%	to a discrete uniform law, between the values (0, pi/2, pi, 3*pi/2).
% 	Such signal is only 'quasi'-analytic.
%	
%	N     : number of points
%	NCOMP : number of points of each component (default: N/5)
% 	F0    : normalized frequency.              (default: 0.25)
% 	Y     : signal
% 	PM0   : initial phase of each component	   (optional).
%
%	Example :
%	 [signal,pm0]=anaqpsk(512,64,0.05); clf; figure(gcf);
%  	 subplot(211); plot(real(signal)); subplot(212); plot(pm0);
%
%	See also ANAFSK, ANABPSK, ANAASK.

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

m=ceil(N/Ncomp);

% jumps=round(3*rand(m,1));
% This is a modification proposed by Alpesh Patel, McMaster University
jumps=floor(4*rand(m,1)); jumps(jumps==4)=3;

pm0=pi*kron(jumps,ones(Ncomp,1))/2; pm0=pm0(1:N,1);
tm=(1:N)'-1;
pm=(2.0*pi*f0*tm+pm0);

y = exp(j*pm);
