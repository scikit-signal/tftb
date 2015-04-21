function S=scale(X,a,fmin,fmax,N,trace);
%SCALE	Scale a signal using the Mellin transform.
%	S=SCALE(X,A,FMIN,FMAX,N,TRACE) computes the A-scaled version
%	of signal X : A^(-1/2) X(T/A) using its Mellin transform.
%
%	X : signal in time to be scaled (Nx=length(X)).
%	A : scale factor. A < 1 corresponds to a compression in the time
%          domain. A can be a vector.		(default : 2)
%	FMIN,FMAX : respectively lower and upper frequency bounds of 
%	   the analyzed signal. These parameters fix the equivalent 
%	   frequency bandwidth (expressed in Hz). When unspecified, you
%	   have to enter them at the command line from the plot of the
%	   spectrum. FMIN and FMAX must be >0 and <=0.5.
%	N : number of analyzed voices (default : automatically determined).
%	TRACE : if nonzero, the progression of the algorithm is shown
%	                                 	(default : 0).
%	S : the A-scaled version of signal X. Length of S can be larger
%          than length of X if A > 1. If A is a vector of length L, S is 
%	   a matrix with L columns. S has the same energy as X.
%
%	Example :
%	 sig=klauder(128); S=scale(sig,2,.05,.45,128);
%	 subplot(211); plot(sig); subplot(212); plot(real(S(65:192)));

%	P. Goncalves, October 1995 - O. Lemoine, June 1996. 
%	Copyright (c) Rice University - CNRS (France)
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

[Mt,xcol] = size(X);

if (nargin == 1),
 a=2;
elseif (nargin==3),
 disp('FMIN will not be taken into account. Determine it with FMAX');
 disp('     from the following plot of the spectrum.'); 
elseif nargin==4,
 N=[];
end;

if nargin <=5,
 trace=0;
end
if (xcol==0)|(xcol>2),
 error('X must have one or two columns');
end; 

Z = hilbert(real(X));
T = Mt;
M = (Mt+rem(Mt,2))/2;

if nargin<=3		        % fmin,fmax,N unspecified
 STF = fft(fftshift(Z)); Nstf=length(STF);
 sp = (abs(STF(1:Nstf/2))).^2; Maxsp=max(sp);
 f = linspace(0,0.5,Nstf/2+1) ; f = f(1:Nstf/2);
 plot(f,sp) ; grid;
 xlabel('Normalized frequency');
 title('Analyzed signal energy spectrum');
 indmin=min(find(sp>Maxsp/1000));
 indmax=max(find(sp>Maxsp/1000));
 fmindflt=max([0.01 0.05*fix(f(indmin)/0.05)]);
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

B = fmax-fmin ; 
R = B/((fmin+fmax)/2) ; 

Nq= ceil((B*T*(1+2/R)*log((1+R/2)/(1-R/2)))/2);
Nmin = Nq-rem(Nq,2);
Ndflt = 2^nextpow2(Nmin);
if nargin<=3,
 Ntxt=['Number of frequency samples (>=',num2str(Nmin),') [',num2str(Ndflt),'] : '];
 N = input(Ntxt);
end
if N~=[],
 if (N<Nmin),
  dispstr=['Warning : the number of analyzed voices (N) should be > ',num2str(Nmin)];
  disp(dispstr);
 end
else
 N=Ndflt; 
end


% Geometric sampling of the analyzed spectrum
k = 1:N;
q = (fmax/fmin)^(1/(N-1));
geo_f  = fmin*(exp((k-1).*log(q)));
t = (1:Mt)-M-1;
tfmatx = zeros(Mt,N);
tfmatx = exp(-2*i*t'*geo_f*pi);
ZS = Z.'*tfmatx; 
ZS(N+1:2*N) = zeros(1,N);


% Mellin transform computation of the analyzed signal
p    = 0:(2*N-1);
MS   = fftshift(ifft(ZS));
beta = (p/N-1)./(2*log(q));


% Inverse Mellin transform and inverse Fourier transform
Mmax = max(ceil(Mt/2*a));
S    = zeros(2*Mmax,length(a)); 
ptr  = 1;
for acurrent = a,
 if trace, disprog(ptr,length(a),10); end
 DMS = exp(-2*i*pi*beta*log(acurrent)).*MS;	% Scaling in Mellin domain
 DS  = fft(fftshift(DMS));			% Inverse Mellin transform
 Mcurrent = ceil(acurrent*Mt/2);
 t = [-Mcurrent:Mcurrent-1]-1;
 itfmatx    = zeros(2*Mcurrent,N);
 itfmatx    = exp(2*i*t'*geo_f*pi);		
 dilate_sig = zeros(2*Mcurrent,1);
 for kk=1:2*Mcurrent,
  dilate_sig(kk) = integ(itfmatx(kk,:).*DS(1:N),geo_f) ;
 end;
 S(Mmax-Mcurrent+1:Mmax+Mcurrent,ptr) = dilate_sig;
 ptr=ptr+1;
end

S=S*norm(X)/norm(S);
