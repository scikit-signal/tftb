function noise=noisecg(N,a1,a2)
%NOISECG Analytic complex gaussian noise.
%	NOISE=NOISECG(N,A1,A2) computes an analytic complex gaussian
%	noise of length N with mean 0.0 and variance 1.0. 
%
%	NOISE=NOISECG(N) yields a complex white gaussian noise.
%
%	NOISE=NOISECG(N,A1) yields a complex colored gaussian noise
%	obtained by filtering a white gaussian noise through a
%		sqrt(1-A1^2)/(1-A1*z^(-1)) 
%	first order filter.
%
%	NOISE=NOISECG(N,A1,A2) yields a complex colored gaussian noise
%	obtained by filtering a white gaussian noise through a
%		sqrt(1-A1^2-A2^2)/(1-A1*z^(-1)-A2*z^(-2)) 
%	second order filter.
%
%	Example :
%	 N=512;noise=noisecg(N);mean(noise),std(noise).^2
%	 subplot(211); plot(real(noise)); axis([1 N -3 3]);
%	 subplot(212); f=linspace(-0.5,0.5,N); 
%	 plot(f,abs(fftshift(fft(noise))).^2);
%		
%	See also RAND, RANDN, NOISECU.

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

if (nargin==1),
 if N<=2,
  noise=(randn(N,1)+j*randn(N,1))/sqrt(2); 
 else
  noise=randn(2^nextpow2(N),1); 
 end
elseif (nargin==2),
 if (abs(a1)>=1.0),
  error ('for a first order filter, abs(a1) must be strictly lower than 1');
 elseif (abs(a1)<=eps),
  if N<=2,
   noise=(randn(N,1)+j*randn(N,1))/sqrt(2); 
  else
   noise=randn(N,1);
  end
 else
  if N<=2,
   noise=(randn(N,1)+j*randn(N,1))/sqrt(2); 
  else
   Nnoise=ceil(N-2.0/log(a1));
   noise=randn(2^nextpow2(Nnoise),1);
  end
  noise=filter(sqrt(1.0-a1^2), [1 -a1],noise);
 end;
elseif (nargin==3),
 if any(roots([1 -a1 -a2])>1),
  error('unstable filter');
 else
  if N<=2,
   noise=(randn(N,1)+j*randn(N,1))/sqrt(2); 
  else
   Nnoise=ceil(N-2.0/log(max(roots([1 -a1 -a2]))));
   noise=randn(2^nextpow2(Nnoise),1);
  end
  noise=filter(sqrt(1.0-a1^2-a2^2), [1 -a1 -a2],noise);
 end;
end;

if N>2,
 noise=hilbert(noise)/std(noise)/sqrt(2);
 noise=noise(length(noise)-(N-1:-1:0));
end

