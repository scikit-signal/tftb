function [tfr,dgr,gam]=tfrgabor(sig,N,q,h,trace)
%TFRGABOR Gabor representation of a signal.
%	[TFR,DGR,GAM]=TFRGABOR(SIG,N,Q,H,TRACE) computes the Gabor
%	representation  of signal X, for a given synthesis window H, on a
%	rectangular grid of size (N,M) in the time-frequency plane. M and N
%	must be such that 
%			N1 = M * N / Q 
%	where N1=length(X) and Q is an integer corresponding to the 
%	degree of oversampling.
%
%	SIG : signal to be analyzed (length(SIG)=N1).
%	N   : number of Gabor coefficients in time (N1 must be a multiple
%	      of N)          	    (default : divider(N1)). 
%	Q   : degree of oversampling ; must be a divider of N 
%				    (default : Q=divider(N)).
%	H   : synthesis window, which was originally chosen as a Gaussian
%	      window by Gabor. Length(H) should be as closed as possible
%	      from N, and must be >=N (default : Gauss(N+1)).
%	      H must be of unit energy, and CENTERED. 
%	TRACE : if nonzero, the progression of the algorithm is shown
%	                                 (default : 0).
%	TFR : Square modulus of the Gabor coefficients. When
%	      called  without output arguments, TFRGABOR runs TFRQVIEW.
%	DGR : Gabor coefficients (complex values). 
%	GAM : biorthogonal (dual frame) window associated to H.
% 
%       If Q=1, the time-frequency plane (TFP) is critically
%	 sampled, so there is no redundancy in the TFP.
%	If Q>1, the TFP is oversampled, allowing a greater
%	 numerical stability of the algorithm. 
%
%	Example :
%	 sig=fmlin(128); 
%	 tfrgabor(sig,64,32); 
%
%	See also all the time-frequency representations listed in
%	 the file CONTENTS (TFR*).

%	O. Lemoine, October 1995 - February 1996.
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
%
%	References : 
%	  o Zibulski, Zeevi "Oversampling in the Gabor Scheme"
%	    IEEE Trans. on Signal Processing, Vol 41, No 8, August 1993
%	    PP 2679-87
%	  o Wexler, Raz "Discrete Gabor Expansions"
%	    Signal Processing, Vol 21, No 3, pp. 207-221, Nov 1990		

N1=length(sig);

if N1<=2,
  error('sig must have at least 3 elements');
elseif nargin<1,
  error('The number of parameters must be at least one');
elseif nargin==1,
  N=divider(N1); q=divider(N); h=tftb_window(odd(N),'Gauss'); trace=0;
elseif nargin==2,
  q=divider(N); h=tftb_window(odd(N),'Gauss'); trace=0;
elseif nargin==3,
  h=tftb_window(odd(N),'Gauss'); trace=0;
elseif nargin==4,
  trace=0;
end
if N<=1|N>=N1,
  error('N must be between 2 and N1-1');
elseif rem(N1,N)~=0,
  error('N must be a divider of N1');
elseif rem(N,q)~=0,
  error('Q must be a divider of N');
end
h=h/norm(h);

% Initializations
if trace, disp('Gabor representation'); end;
[xrow,xcol]=size(sig);
[wrow,wcol]=size(h);
if xrow>xcol,
 sig=sig.';
elseif wrow>wcol,
 h=h.';
end

M=q*N1/N;
Mb=N1/N;        % frequency-shift between two Gabor coefficients
Nb=N1/M;        % time-shift between two Gabor coefficients

% Zak transform of h
Nh=length(h);
if rem(Nh,2)==0,
 error('H must have an odd length');
end

% Time shift to be correctly localized in time :
alpha=round((2*N1/N-1-Nh)/(2*q));
hN1=zeros(N1,1); 
hN1(round((N1-(Nh-1))/2)-alpha:round((N1+Nh-1)/2)-alpha)=h;	

DZTh=zeros(M,Nb);
Msig=reshape(hN1,Nb,M);
DZTh=fft(Msig.')/sqrt(M);

% Sum of elements of h-zak transform

Mzh=zeros(M,Nb);
x=(1:M);

for l=0:q-1,
 mod=modulo(x-l*M/q,M);
 Mzh=Mzh+abs(DZTh(mod,:)).^2;
end

[ind1,ind2]=find(Mzh<eps);
for l=1:length(ind1),
 Mzh(ind1(l),ind2(l))=1;
end
 
% Zak transform of the biorthogonal (dual frame) window gam

DZTgam=zeros(M,Nb);
DZTgam=DZTh./Mzh;
gam=real(izak(DZTgam))/N1;

% Computation of the Gabor coefficients from the dual frame window

dgr=zeros(M,N);
tfr=zeros(M,N);
k=1:N1;
dgrN1=zeros(N1,N);
 
for n=1:N, 
 if trace, disprog(n,N,10); end;
 indice=modulo(k-n*M/q,N1);
 dgrN1(:,n)=fft(sig.*fftshift(gam(indice)')).';
end

dgr=dgrN1(1:Nb:N1,:);
tfr=abs(dgr).^2;
t=1:M/q:N1;

if (nargout==0),
 tfrqview(tfr,sig.',t,'tfrgabor',N,q,h);
end;
