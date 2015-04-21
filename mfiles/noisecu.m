function noise=noisecu(N);
%NOISECU Analytic complex uniform white noise.
%	NOISE=NOISECU(N) computes an analytic complex white uniform
%	 noise of length N with mean 0.0 and variance 1.0. 
%
%	Examples :
%	 N=512;noise=noisecu(N);mean(noise),std(noise).^2
%	 subplot(211); plot(real(noise)); axis([1 N -1.5 1.5]);
%	 subplot(212); f=linspace(-0.5,0.5,N); 
%	 plot(f,abs(fftshift(fft(noise))).^2);
%		
%	See also RAND, RANDN, NOISECG.

%	O. Lemoine, June 95/May 96 - F. Auger, August 95.
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


if (N <= 0),
 error ('The signal length N must be strictly positive' );
end;

MatlabVersion=version; MatlabVersion=str2num(MatlabVersion(1));
% I hope that future Matlab versions will be more compatible
if (MatlabVersion==4),
 rand('uniform');
elseif (MatlabVersion>5),
 error('unsupported matlab version. please send an email.'); 
end;


if N<=2,
 noise=(rand(N,1)-0.5+j*(rand(N,1)-0.5))*sqrt(6); 
else
 noise=rand(2^nextpow2(N),1)-0.5;
end

if N>2,
 noise=hilbert(noise)/std(noise)/sqrt(2);
 noise=noise(length(noise)-(N-1:-1:0));
end
