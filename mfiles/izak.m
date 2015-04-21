function sig=izak(dzt)
%IZAK 	Inverse Zak transform.
%	SIG=IZAK(DZT) computes the inverse Zak transform of matrix DZT.
%
%	DZT : (N,M) matrix of Zak samples.
%	SIG : Output signal (M*N,1) containing the inverse Zak transform.
%
%	Example :
%	 sig=fmlin(256); DZT=zak(sig); sigr=izak(DZT);
%	 plot(real(sigr)); hold; plot(real(sig)); hold;
%
%	See also ZAK, TFRGABOR.

%	O. Lemoine - February 1996
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

[N,M]=size(dzt);

if nargin<1,
  error('The number of parameters must be one');
end

sig=zeros(N*M,1);

for m=1:M,
  sig(m+(0:N-1)*M)=sqrt(N)*ifft(dzt(:,m));
end

