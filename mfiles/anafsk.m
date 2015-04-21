function [y,iflaw]=anafsk(N,Ncomp,Nbf);
%ANAFSK Frequency Shift Keying (FSK) signal.
%	[Y,IFLAW]=ANAFSK(N,NCOMP,NBF) simulates a phase coherent
%	Frequency Shift Keying (FSK) signal. This signal is a succession
%	of complex sinusoids of NCOMP points each and with a normalized
%	frequency uniformly chosen between NBF distinct values between
%	0.0 and 0.5.  
% 	Such signal is only 'quasi'-analytic.
%	
%	N     : number of points
%	NCOMP : number of points of each component (default: N/5)
%	NBF   : number of distinct frequencies     (default: 4  )
% 	Y     : signal
% 	IFLAW : instantaneous frequency law  (optional).
%
%	Example :
%	 [signal,ifl]=anafsk(512,64,5); clf; figure(gcf);
%  	 subplot(211); plot(real(signal)); subplot(212); plot(ifl);
%
%	See also ANABPSK, ANAQPSK, ANAASK.

%	O. Lemoine - June 1995, F. Auger - August 1995.
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
 Ncomp=round(N/5); Nbf=4;
elseif (nargin == 2),
 Nbf=4;
end;

if (N <= 0),
 error ('The signal length N must be strictly positive' );
end;

MatlabVersion=version; MatlabVersion=str2num(MatlabVersion(1));
if (MatlabVersion==4), rand('uniform'); end;

m=ceil(N/Ncomp);
freqs=0.25+0.25*(floor(Nbf*rand(m,1))/Nbf-(Nbf-1)/(2*Nbf));
iflaw=kron(freqs,ones(Ncomp,1)); iflaw=iflaw(1:N,1);
y=exp(j*2*pi*cumsum(iflaw));

