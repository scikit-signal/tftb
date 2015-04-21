function [fm,B]=locfreq(sig);
%LOCFREQ Frequency localization caracteristics.
%	[FM,B]=LOCFREQ(SIG) computes the frequency localization 
%	caracteristics of signal SIG.
% 
%	SIG   is the signal.
%	FM    is the averaged normalized frequency center.
%	B     is the frequency spreading.
%
%	Example :
%	 z=amgauss(160,80,50);[tm,T]=loctime(z),[fm,B]=locfreq(z),B*T
%
%	See also LOCTIME.

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

[N,sigc]=size(sig);
if (sigc~=1),
 error('The signal must have 1 column');
else
 No2r=round(N/2);
 No2f=fix(N/2);
 Sig=fft(sig);
 Sig2=abs(Sig).^2;
 Sig2=Sig2/mean(Sig2);
 freqs=[0:No2f-1 -No2r:-1]'/N;
 fm=mean(freqs.*Sig2);
 B=2*sqrt(pi*mean((freqs-fm).^2.*Sig2));
end;
