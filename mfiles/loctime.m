function [tm,T]=loctime(sig);
%LOCTIME Time localization caracteristics.
%	[TM,T]=LOCTIME(SIG) computes the time localization
%	caracteristics of signal SIG. 
% 
%	SIG is the signal.
%	TM  is the averaged time center.
%	T   is the time spreading.
%
%	Example :
%	 z=amgauss(160,80,50); [tm,T]=loctime(z)
%
%	See also LOCFREQ.

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

[sigr,sigc]=size(sig);
if (sigc~=1),
 error('The signal must have 1 column');
else
 sig2=abs(sig).^2; sig2=sig2/mean(sig2);
 t=(1:sigr)';
 tm=mean(t.*sig2);
 T=2*sqrt(pi*mean((t-tm).^2 .* sig2)); 
end;
