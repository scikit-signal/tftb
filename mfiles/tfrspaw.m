function [tfr,t,f]=tfrspaw(X,time,K,nh0,ng0,fmin,fmax,N,trace);
%TFRSPAW Smoothed Pseudo Affine Wigner time-frequency distributions.
%	[TFR,T,F]=TFRSPAW(X,T,K,NH0,NG0,FMIN,FMAX,N,TRACE)
%	generates the auto- or cross- Smoothed Pseudo Affine Wigner
%	distributions.  
%
%	X : signal (in time) to be analyzed. If X=[X1 X2], TFRSPAW 
%	   computes the cross-Smoothed Pseudo Affine Wigner distribution.
%						(Nx=length(X)).
%	T : time instant(s) on which the TFR is evaluated  
%	   					(default : 1:Nx).
%	K : label of the K-Bertrand distribution. The distribution with
%	   parameterization function 
%	   lambdak(u,K) = (K (exp(-u)-1)/(exp(-Ku)-1))^(1/(K-1)) 
%	   is computed				(default : 0).
%	     K=-1 : Smoothed pseudo (active) Unterberger distribution 
%	     K=0  : Smoothed pseudo Bertrand distribution
%	     K=1/2: Smoothed pseudo D-Flandrin distribution
%	     K=2  : Affine smoothed pseudo Wigner-Ville distribution.
%	NH0 : half length of the analyzing wavelet at coarsest scale.  
%	   A Morlet wavelet is used. NH0 controles the frequency 
%	   smoothing of the smoothed pseudo Affine Wigner distribution.
%						(default : sqrt(Nx)).
%	NG0 : half length of the time smoothing window. 
%	   NG0 = 0 corresponds to the Pseudo Affine Wigner distribution.  
%						(default : 0).
%	FMIN,FMAX : respectively lower and upper frequency bounds of 
%	   the analyzed signal. These parameters fix the equivalent 
%	   frequency bandwidth (expressed in Hz). When unspecified, you
%	   have to enter them at the command line from the plot of the
%	   spectrum. FMIN and FMAX must be >0 and <=0.5. 
%	N : number of analyzed voices (default : automatically determined).
%	TRACE : if nonzero, the progression of the algorithm is shown
%						(default : 0).
%	TFR : time-frequency matrix containing the coefficients of the
%	   decomposition (abscissa correspond to uniformly sampled time,
%	   and ordinates correspond to a geometrically sampled
%	   frequency). First row of TFR corresponds to the lowest 
%	   frequency. When called without output arguments, TFRSPAW
%	   runs TFRQVIEW.
%	F : vector of normalized frequencies (geometrically sampled 
%	   from FMIN to FMAX).
%
%	Example :    
%	 sig=altes(64,0.1,0.45); tfrspaw(sig);
%	 
%	See also all the time-frequency representations listed in
%	the file CONTENTS (TFR*)

%	P. Goncalves, October 95 - O. Lemoine, June 1996.
%	Copyright (c) 1995 Rice University - CNRS (France) 1996.
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
if (nargin<=8), trace=0; end

if (nargin == 1),
 time=1:xrow; K=0; nh0=sqrt(xrow); ng0=0;
elseif (nargin == 2),
 K=0; nh0=sqrt(xrow); ng0=0;
elseif (nargin == 3),
 nh0=sqrt(xrow); ng0=0;
elseif (nargin == 4),
 ng0=0;
elseif (nargin == 6),
 disp('FMIN will not be taken into account. Determine it with FMAX');
 disp('     from the following plot of the spectrum.'); 
elseif (nargin == 7),
 N=[];
end;

[trow,tcol] = size(time);
if (xcol==0)|(xcol>2),
 error('X must have one or two columns');
elseif (trow~=1),
 error('T must only have one row'); 
end; 

Mt=length(X); 

if trace, 
 if (K==-1),
  disp('Smoothed pseudo (active) Unterberger distribution');
 elseif K==0,
   disp('Smoothed pseudo Bertrand distribution');
 elseif K==1/2,
   disp('Smoothed pseudo D-Flandrin distribution');
 elseif K==2,
   disp('Affine smoothed pseudo Wigner-Ville distribution');
 else
  disp('Smoothed Pseudo Affine Wigner distribution');
 end;
end;

if xcol==1,
 X1=X;
 X2=X; 
else
 X1=X(:,1);
 X2=X(:,2);
end
s1 = real(X1);
s2 = real(X2);
M  = (Mt+rem(Mt,2))/2;

if nargin<=6,				        % fmin,fmax,N unspecified
 STF1 = fft(fftshift(s1(min(time):max(time))));
 Nstf=length(STF1);
 sp1 = (abs(STF1(1:Nstf/2))).^2; Maxsp1=max(sp1);
 STF2 = fft(fftshift(s2(min(time):max(time)))); 
 sp2 = (abs(STF2(1:Nstf/2))).^2; Maxsp2=max(sp2);
 f = linspace(0,0.5,Nstf/2+1) ; f=f(1:Nstf/2);
 plot(f,sp1) ; grid; hold on ; plot(f,sp2) ; hold off
 xlabel('Normalized frequency');
 title('Analyzed signal energy spectrum');
 axis([0 1/2 0 1.2*max(Maxsp1,Maxsp2)]) ; 
 indmin=min(find(sp1>Maxsp1/100));
 indmax=max(find(sp1>Maxsp1/100));
 fmindflt=max([0.01 0.05*fix(f(indmin)/0.05)]);
 fmaxdflt=0.05*ceil(f(indmax)/0.05);
 txtmin=['Lower frequency bound [',num2str(fmindflt),'] : '];
 txtmax=['Upper frequency bound [',num2str(fmaxdflt),'] : '];
 fmin = input(txtmin); fmax = input(txtmax);
 if isempty(fmin), fmin=fmindflt; end
 if isempty(fmax), fmax=fmaxdflt; end
end

if (fmin >= fmax)
 error('FMAX must be greater or equal to FMIN');
elseif fmin<=0.0 | fmin>0.5,
 error('FMIN must be > 0 and <= 0.5');
elseif fmax<=0.0 | fmax>0.5,
 error('FMAX must be > 0 and <= 0.5');
end

B    = fmax-fmin ; 
R    = B/((fmin+fmax)/2) ; 
Qte  = fmax/fmin ;    
umax = log(Qte); 
Teq  = nh0/(fmax*umax);  
if Teq<2*nh0,
 M0 = (2*nh0^2)/Teq-nh0+1;
else
 M0 = 0;
end;
MU = round(nh0+M0);
T  = 2*MU-1;

Nq = ceil((B*T*(1+2/R)*log((1+R/2)/(1-R/2)))/2);
Nmin = Nq-rem(Nq,2);
Ndflt = 2^nextpow2(Nmin);
if nargin<=6,
 Ntxt=['Number of frequency samples (>=',num2str(Nmin),') [',num2str(Ndflt),'] : '];
 N = input(Ntxt);
end
if ~isempty(N),
 if (N<Nmin),
  dispstr=['Warning : the number of analyzed voices (N) should be > ',num2str(Nmin)];
  disp(dispstr);
 end
else
 N=Ndflt; 
end

fmin_s = num2str(fmin); 
fmax_s = num2str(fmax); 
N_s = num2str(N);
if trace,
 disp(['Frequency runs from ',fmin_s,' to ',fmax_s,' with ',N_s,' points']);
end

k = 1:N;
q = (fmax/fmin)^(1/(N-1));
a = exp((k-1).*log(q));         % a is an increasing scale vector.
geo_f = fmin*a;                 % geo_f is a geometrical increasing
                                % frequency vector.

% Wavelet decomposition computation
matxte1 = zeros(N,tcol);
matxte2 = zeros(N,tcol);
[p1,p2,p3,wt1] = tfrscalo(s1,time,nh0,fmin,fmax,N) ;
[p1,p2,p3,wt2] = tfrscalo(s2,time,nh0,fmin,fmax,N) ;
for ptr = 1:N, 
 matxte1(ptr,:) = wt1(ptr,:).*sqrt(a(N-ptr+1)) ; 
 matxte2(ptr,:) = wt2(ptr,:).*sqrt(a(N-ptr+1)) ; 
end ;

umin = -umax;
u=linspace(umin,umax,2*MU+1);
du = u(2)-u(1);
u=u(1:2*MU);
u(MU+1) = 0;
p = 0:(2*N-1);
beta = (p/N-1)./(2*log(q));
l1=zeros(2*MU,2*N);
l2=zeros(2*MU,2*N);
for m = 1:2*MU,
 l1(m,:) = exp(-2*i*pi*beta*log(lambdak( u(m),K)));
 l2(m,:) = exp(-2*i*pi*beta*log(lambdak(-u(m),K)));
end 


% Determination of the time smoothing window G
if ng0==0,
 G = ones(2*MU,1);
else
 a_t = 3 ;            % (attenuation of 10^(-a_t) at t = tmax)
 sigma_t = ng0*fmax/sqrt(2*a_t*log(10));
 a_u = 2 * pi^2 * sigma_t^2 * umax^2 / log(10) ;
 sigma_u = 1/(2 * pi * sigma_t) ;
 G = exp(-(a_u*log(10)/MU^2)*[-MU:MU-1].^2); 
 if sigma_u < du
  disp('Maximum time smoothing reached. Increase width of wavelet for effectiveness.') ;
 end
 G=G';
end

waf = zeros(2*MU,N);
tfr = zeros(N,tcol);
S1  = zeros(1,2*N);
S2  = zeros(1,2*N);
MX1 = zeros(2*N,2*MU);
MX2 = zeros(2*N,2*MU);
TX1 = zeros(2*MU,N);
TX2 = zeros(2*MU,N);

for ti = 1:tcol,
 if trace, disprog(ti,tcol,10); end

 S1(1:N) = matxte1(:,ti).';
 Mellin1 = fftshift(ifft(S1));
 MX1 = (l1.*Mellin1(ones(1,2*MU),:)).';
 MX1 = fft(MX1);
 TX1 = MX1(1:N,:).';

 S2(1:N) = matxte2(:,ti).';
 Mellin2 = fftshift(ifft(S2));     
 MX2 = (l2.*Mellin2(ones(1,2*MU),:)).';
 MX2 = fft(MX2);
 TX2 = MX2(1:N,:).';

 waf = real(TX1.*conj(TX2)).*G(:,ones(N,1));
 tfr(:,ti) = (sum(waf).*geo_f).';	% first row of tfr corresponds to
				% the lowest frequency.
end;


t = time; 
f = geo_f'; 

% Normalization
SP1 = fft(hilbert(s1)); 
SP2 = fft(hilbert(s2)); 
indmin = 1+round(fmin*(xrow-2));
indmax = 1+round(fmax*(xrow-2));
SP1ana = SP1(indmin:indmax);
SP2ana = SP2(indmin:indmax);
tfr = tfr*(SP1ana'*SP2ana)/integ2d(tfr,t,f)/N;


if (nargout==0),
 tfrqview(real(tfr),hilbert(real(X)),t,'tfrspaw',K,nh0,ng0,N,f);
end;

