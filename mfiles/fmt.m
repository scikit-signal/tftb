function [mellin,beta]=fmt(X,fmin,fmax,N);
%FMT    Fast Fourier Mellin Transform.
%       [MELLIN,BETA]=FMT(X,FMIN,FMAX,N) computes the Fast Mellin
%       Transform of signal X.
%
%       X : signal in time (Nx=length(X)).
%       FMIN,FMAX : respectively lower and upper frequency bounds of 
%        the analyzed signal. These parameters fix the equivalent 
%        frequency bandwidth (expressed in Hz). When unspecified, you
%        have to enter them at the command line from the plot of the
%        spectrum. FMIN and FMAX must be >0 and <=0.5.         
%       N : number of analyzed voices. N must be even
%				 (default : automatically determined).
%       MELLIN : the N-points Mellin transform of signal S.
%       BETA : the N-points Mellin variable.
%
%       Examples :   
%        sig=altes(128,0.05,0.45); 
%	 [MELLIN,BETA]=fmt(sig,0.05,0.5,128);
%        plot(BETA,real(MELLIN));
%
%       See also IFMT, FFT, IFFT.

%       P. Goncalves 9-95 - O. Lemoine, June 1996. 
%       Copyright (c) 1995 Rice University - CNRS (France) 1996.
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
 error('At least one parameter required');
end;

[xrow,xcol] = size(X);

if (nargin==2),
 disp('FMIN will not be taken into account. Determine it with FMAX');
 disp('     from the following plot of the spectrum.'); 
elseif nargin==3,
 N=[];
elseif (nargin==4 & rem(N,2)~=0),
 error('N must be even');
end;
if (xcol==0)|(xcol>2),
 error('X must have one or two columns');
end

Mt = length(X); 
Z  = hilbert(real(X));
M  = (Mt+rem(Mt,2))/2;

if nargin<=2,                                   % fmin,fmax,N unspecified
 STF = fft(fftshift(X)); Nstf=length(STF);
 sp = (abs(STF(1:Nstf/2))).^2; Maxsp=max(sp);
 f = linspace(0,0.5,Nstf/2+1) ; f = f(1:Nstf/2);
 plot(f,sp) ; grid;
 xlabel('Normalized frequency');
 title('Analyzed signal energy spectrum');
 indmin=min(find(sp>Maxsp/1000));
 indmax=max(find(sp>Maxsp/1000));
 fmindflt=max([0.001 0.05*fix(f(indmin)/0.05)]);
 fmaxdflt=0.05*ceil(f(indmax)/0.05);
 txtmin=['Lower frequency bound [',num2str(fmindflt),'] : '];
 txtmax=['Upper frequency bound [',num2str(fmaxdflt),'] : '];
 fmin = input(txtmin); fmax = input(txtmax);
 if fmin==[], fmin=fmindflt; end
 if fmax==[], fmax=fmaxdflt; end
end

if fmin >= fmax
 error('FMAX must be greater or equal to FMIN');
elseif fmin<=0.0 | fmin>0.5,
 error('FMIN must be > 0 and <= 0.5');
elseif fmax<=0.0 | fmax>0.5,
 error('FMAX must be > 0 and <= 0.5');
end

B = fmax-fmin;       		% Bandwidth of the signal X
R = B/((fmin+fmax)/2);		% Relative bandwidth of X

Nq= ceil((B*Mt*(1+2/R)*log((1+R/2)/(1-R/2)))/2);
Nmin = Nq-rem(Nq,2);
Ndflt = 2^nextpow2(Nmin);
if nargin<=2,
 Ntxt=['Number of frequency samples (>=',num2str(Nmin),') [',num2str(Ndflt),'] : '];
 N = input(Ntxt);
end

if N~=[],
 if (N<Nmin),
  dispstr=['Warning : the number of analyzed voices (N) should be >= ',num2str(Nmin)];
  disp(dispstr);
 end
else
 N=Ndflt; 
end


% Geometric sampling of the analyzed spectrum
No2 = N/2;
k = (1:No2);
q = (fmax/fmin)^(1/(No2-1));
t = (1:Mt)-M-1;
geo_f  = fmin*(exp((k-1).*log(q)));
tfmatx = zeros(Mt,N);
tfmatx = exp(-2*j*pi*t'*geo_f);
ZS = Z.'*tfmatx; 
ZS(No2+1:N) = zeros(1,N-No2);


% Mellin transform computation of the analyzed signal
p = 0:(N-1);
L = log(fmin)/log(q);
mellin = N*log(q)*fftshift(ifft(ZS)).*exp(j*2*pi*L*(p/N-1/2));
beta   = (p/N-1/2)./log(q);


% Normalization
SP = fft(hilbert(real(X))); 
indmin = 1+round(fmin*(xrow-2));
indmax = 1+round(fmax*(xrow-2));
SPana = SP(indmin:indmax);
nu = (indmin:indmax)'/N; 
SPp = SPana./nu;
Normsig = sqrt(SPp'*SPana);

mellin = mellin*Normsig/norm(mellin);
