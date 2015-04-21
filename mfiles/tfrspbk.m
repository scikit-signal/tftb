function [tfr,t,f]=tfrspbk(X,time,K,nh0,ng0,fmin,fmax,N,trace);
%TFRSPBK Smoothed Pseudo K-Bertrand time-frequency distribution.
%	[TFR,T,F]=TFRSPBK(X,T,K,NH0,NG0,FMIN,FMAX,N,TRACE)
%	generates the auto- or cross- Smoothed Pseudo K-Bertrand
%	distribution.  
%
%	X : signal (in time) to be analyzed. If X=[X1 X2], TFRSPBK 
%	   computes the cross-Smoothed Pseudo K-Bertrand distribution.
%						(Nx=length(X)).
%	T : time instant(s) on which the TFR is evaluated. TIME must 
%	   be a uniformly sampled vector whose elements are between 1 
%	   and Nx.				(default : 1:Nx).
%	K : label of the K-Bertrand distribution. The distribution with
%	   parametrization function 
%	   lambdak(u,K) = (K (exp(-u)-1)/(exp(-Ku)-1))^(1/(K-1)) 
%	   is computed				(default : 0).
%	     K=-1 : Smoothed pseudo (active) Unterberger distribution 
%	     K=0  : Smoothed pseudo Bertrand distribution
%	     K=1/2: Smoothed pseudo D-Flandrin distribution
%	     K=2  : Affine smoothed pseudo Wigner-Ville distribution.
%	NH0 : half length of the analyzing wavelet at coarsest scale.  
%	   A Morlet wavelet is used. NH0 controles the frequency 
%	   smoothing of the smoothed pseudo K-Bertrand distribution.
%						(default : sqrt(Nx)).
%	NG0 : half length of the time smoothing window. 
%	   NG0 = 0 corresponds to the Pseudo K-Bertrand distribution.  
%						(default : 0).
%	FMIN,FMAX : respectively lower and upper frequency bounds of 
%	   the analyzed signal. These parameters fix the equivalent 
%	   frequency bandwidth (expressed in Hz). When unspecified, you
%	   have to enter them at the command line from the plot of the
%	   spectrum. FMIN and FMAX must be >0 and <=0.5. 
%	N : number of analyzed voices	 	(default : Nx).
%	TRACE : if nonzero, the progression of the algorithm is shown
%						(default : 0).
%	TFR : time-frequency matrix containing the coefficients of the
%	   decomposition (abscissa correspond to uniformly sampled time,
%	   and ordinates correspond to a geometrically sampled
%	   frequency). First row of TFR corresponds to the lowest 
%	   frequency. When called without output arguments, TFRSPBK
%	   runs TFRQVIEW.
%	F : vector of normalized frequencies (geometrically sampled 
%	    from FMIN to FMAX).
%
%	Example :    
%	 sig=altes(64,0.1,0.45); tfrspbk(sig);
%	 
%	See also TFRBERT, TFRUNTAC, TFRUNTPA, TFRSCALO, TFRDFLA, TFRASPW.

%	P. Goncalves, October 95 - O. Lemoine, June 1996.
%	Copyright (c) 1995 Rice University
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
if nargin<=8, trace=0; end

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
 N=xrow;
end;

[trow,tcol] = size(time);
if (xcol==0)|(xcol>2),
 error('X must have one or two columns');
elseif (trow~=1),
 error('TIME must only have one row'); 
end; 

Mt=length(X); 

if trace, disp('Smoothed Pseudo K-Bertrand distribution'); end;

if xcol==1,
 X1=X;
 X2=X; 
else
 X1=X(:,1);
 X2=X(:,2);
end
s1 = real(X1);
s2 = real(X2);

if rem(Mt,2)~=0, 
 s1 = [s1;0];
 s2 = [s2;0]; 
 M  = (Mt+1)/2;
else
 M  = Mt/2;
end ;

t = [-nh0:nh0-1];
Tmin = 1 ;
Tmax = 2*nh0 ; 
T = Tmax-Tmin ;

if nargin<=6,				        % fmin,fmax,N unspecified
 STF1 = fft(fftshift(s1(min(time):max(time)))); Nstf=length(STF1);
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
Qte = fmax/fmin ;    
umax = log(Qte); 
Teq = nh0/(fmax*umax);  
if Teq<2*nh0,
 M0 = round((2*nh0^2)/Teq-nh0)+1;
 MU = nh0+M0;
 T2 = 2*MU-1;
else
 M0 = 0;
 MU = nh0;
 T2 = 2*MU-1;
end;

if nargin<=6,
 Nq= ceil((B*T2*(1+2/R)*log((1+R/2)/(1-R/2)))/2);
 Nmin = Nq-rem(Nq,2);
 Ndflt = 2^nextpow2(Nmin);
 Ntxt=['Number of frequency samples (>=',num2str(Nmin),') [',num2str(Ndflt),'] : '];
 N = input(Ntxt);
 if N==[], N=Ndflt; end
end

fmin_s = num2str(fmin) ; fmax_s = num2str(fmax) ; N_s = num2str(N) ;
if trace,
 disp(['Frequency runs from ',fmin_s,' to ',fmax_s,' with ',N_s,' points']);
end


k = 1:N;
q = (fmax/fmin)^(1/(N));
a = (exp((k-1).*log(q)));       % a is an increasing scale vector.
geo_f(k) = fmin*a ;             % geo_f is a geometrical increasing
                                % frequency vector.

% Morlet wavelet decomposition computation
t0 = 1; 
t1 = Mt; 
Mtr = Mt;
z1 = hilbert(s1).';
matxte1 = zeros(N,Mt);
z2 = hilbert(s2).';
matxte2 = zeros(N,Mt);
nu0 = geo_f(N);
for ptr=1:N,
 nha = round(nh0*a(ptr));
 nua = nu0/a(ptr);
 ha = exp(-(2*log(10)/(nh0*a(ptr))^2)*(-nha:nha).^2).*exp(-i*2*pi*nua*(-nha:nha));
 detail1 = conv(z1(t0:t1),fliplr(ha)); 
 matxte1(N-ptr+1,:) = detail1(nha+1:length(detail1)-nha);
 detail2 = conv(z2(t0:t1),fliplr(ha)); 
 matxte2(N-ptr+1,:) = detail2(nha+1:length(detail2)-nha);
            % first row of matxte corresponds to the lowest frequency.
end;


% Pseudo-Bertrand distribution computation
tfr=zeros(N,tcol);
umin = -umax;
u=linspace(umin,umax,2*MU+1);
U(MU+1) = 0;
k = 1:2*N;
beta(k) = -1/(2*log(q))+(k-1)./(2*N*log(q));
for m = 1:2*MU+1,
 l1(m,:) = exp(-(2*i*pi*beta+1/2).*log(lambdak(u(m),K)));
end 
if ng0==0
 decay = 0 ;
elseif ng0~=0
 gamma0 = ng0*fmax ;
 alpha = - log(0.01)/gamma0^2 ; 
 u0 = sqrt(-alpha*log(-0.01*sqrt(alpha/pi))/pi^2) ;  
 decay = -log(0.01)*(umax/u0)^2/log(10) ; 
end
if decay==Inf 
 G = zeros(1,2*MU) ;
 G(MU+1) = 1 ;
elseif decay==0,
 G=ones(1,2*MU);
else
 Nb=2*MU+1;
 G = amgauss(Nb,(Nb+1)/2,(Nb-1)*sqrt(pi/(decay*log(10)))/2).';
 G = G(1:2*MU) ;
end
xx = exp(-[0:N-1].*log(q));
xx = xx(ones(1,2*MU),:).*G(ones(1,N),:)';
indi=1;
for ti = time,
 if trace, disprog(ti-time(1)+1,time(tcol)-time(1)+1,10); end
 S1 = zeros(1,2*N);
 S1(1:N) = matxte1(:,ti).';
 Mellin1 = fftshift(ifft(S1.*exp([0:2*N-1].*log(q)))) ;
 S2 = zeros(1,2*N);
 S2(1:N) = matxte2(:,ti).';
 Mellin2 = fftshift(ifft(S2.*exp([0:2*N-1].*log(q)))) ;     
 waf = zeros(2*MU,N) ; 
 MX1 = l1.*Mellin1(ones(1,2*MU+1),:) ;
 X1 = fft(MX1.');
 X1 = X1(1:N,:).' ;
 MX2 = l1.*Mellin2(ones(1,2*MU+1),:) ;
 X2 = fft(MX2.');
 X2 = X2(1:N,:).';      
 waf = real(X1(1:2*MU,:).*conj(X2(2*MU+1:-1:2,:)).*xx) ;
 tfr(:,indi) = sum(waf).';		% first row of tfr corresponds to
 indi = indi+1; 			% the lowest frequency.
end;

t = time; 
f = geo_f';
tfr = tfr./integ2d(tfr,t,f)*sum(s1.*conj(s2)) ;

disp(' ');

if (nargout==0),
 tfrqview(tfr,X,t,'tfrspbk',K,nh0,ng0,N,f);
end;




