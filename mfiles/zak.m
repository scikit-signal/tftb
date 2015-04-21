function dzt=zak(sig,N,M)
%ZAK 	Zak transform.
%	DZT=ZAK(SIG,N,M) computes the Zak transform of signal SIG.
%
%	SIG : Signal to be analysed (length(X)=N1).
%	N   : number of Zak coefficients in time (N1 must be a multiple
%	      of N)          (default : closest integer towards sqrt(N1)). 
%	M   : number of Zak coefficients in frequency (N1 must be a
%	      multiple of M) (default : M=N1/N).
%	DZT : Output matrix (N,M) containing the discrete Zak transform.
%
%	Example : 
%	 sig=fmlin(256); DZT=zak(sig);
%	 imagesc(DZT);
%
%	See also IZAK, TFRGABOR.

%	O. Lemoine - February 1996.
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

N1=length(sig);

if N1<=2,
  error('SIG must have at least 3 elements');
elseif nargin<1,
  error('The number of parameters must be at least one');
elseif nargin==1,
  [N,M]=divider(N1);
elseif nargin==2,
  N=N1/M;
elseif M<=1|N<=1|M>=N1|N>=N1,
  error('M and N must be between 2 and N1-1');
elseif rem(N1,M)~=0|rem(N1,N)~=0,
  error('M and N must be such that N1 is a multiple of M and N');
end

Msig=zeros(M,N);
dzt=zeros(N,M);

Msig=reshape(sig,M,N);

dzt=fft(Msig.')/sqrt(N);

   