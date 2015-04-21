function [d,f]=d2statio(sig);
%D2STATIO Distance to stationarity
%	[D,F]=D2STATIO(SIG) evaluates the distance of the signal
%	to stationarity, using the pseudo Wigner-Ville distribution.
%
%	SIG : signal to be analyzed (real or complex).
%	D   : vector giving the distance to stationarity for each frequency.
%	F   : vector of frequency bins
%
%	Example :
%	 sig=noisecg(128); [d,f]=d2statio(sig); plot(f,d);
%	 xlabel('Frequency'); ylabel('Distance'); 
%
%	 sig=fmconst(128); [d,f]=d2statio(sig); plot(f,d);
%	 xlabel('Frequency'); ylabel('Distance'); 
%
 
%	O. Lemoine - May 1996.
%	Copyright (c) by CNRS France, 1996.
%       send bugs to f.auger@ieee.org

N=length(sig);

[tfr,t,f]=tfrspwv(sig);
 
d2=((tfr-mean(tfr')'*ones(1,N))/norm(sig)).^2;	

d=mean(d2')';

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
