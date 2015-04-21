function x=ifmt(mellin,beta,M);
%IFMT	Inverse fast Mellin transform.
%	X=IFMT(MELLIN,BETA,M) computes the inverse fast Mellin
%	transform of MELLIN.
%	WARNING : the inverse of the Mellin transform is correct only 
%	if the Mellin transform has been computed from FMIN to 0.5 Hz, 
%	and if the original signal is analytic.
%	
%	MELLIN : Mellin transform to be inverted. Mellin must have been
%	 obtained from FMT with frequency running from FMIN to 0.5 Hz.
%	BETA : Mellin variable issued from FMT.
%	M : number of points of the inverse Mellin transform. 
%					(default : length(MELLIN)).
%	X : inverse Mellin transform with M points in time.
%
%	Example :   
%	 sig=atoms(128,[64,0.25,32,1]); 
%	 [MELLIN,BETA]=fmt(sig,0.05,0.5,256); clf;
%	 X=ifmt(MELLIN,BETA,128); plot(real(X)); hold; 
%	 plot(real(sig),'g'); hold; 
%
%	See also : fmt, fft, ifft.
%

%	P. Goncalves 9-95 - O. Lemoine, June 1996. 
%	Copyright (c) 1995 Rice University - CNRS (France).
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

N=length(mellin);

if nargin<2,
 error('At least 2 input parameters required');
elseif nargin==2,
 M=N;
end

No2 = (N+rem(N,2))/2;
q   = exp(1/(N*(beta(2)-beta(1))));
fmin = 0.5/(q^(No2-1));


% Inverse Mellin transform computation 
p = 0:(N-1);
L = log(fmin)/log(q);
S = fft(fftshift(mellin.*exp(-j*2*pi*L*(p/N-1/2))/(N*log(q))));
S = S(1:No2);


% Inverse Fourier transform
k = (1:No2);
x = zeros(M,1); 
t = (1:M)-(M+rem(M,2))/2-1;
geo_f = fmin*(exp((k-1).*log(q))) ;
itfmatx = zeros(M,No2);
itfmatx = exp(2*i*t'*geo_f*pi);
for k=1:M,
 x(k) = real(integ(itfmatx(k,:).*S,geo_f));
end;

x = hilbert(x);

% Normalization
SP = fft(x); fmax=0.5;
indmin = 1+round(fmin*(M-2));
indmax = 1+round(fmax*(M-2));
SPana = SP(indmin:indmax);
nu = (indmin:indmax)'/M; 
SPp = SPana./nu;
Esm = SPp'*SPana;
x = x*norm(mellin)/sqrt(Esm);



