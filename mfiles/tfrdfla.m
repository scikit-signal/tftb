function  [tfr,t,f]=tfrdfla(X,time,fmin,fmax,N,trace);
%TFRDFLA D-Flandrin time-frequency distribution.
%	[TFR,T,F]=TFRDFLA(X,T,FMIN,FMAX,N,TRACE) generates the
%	auto- or cross- D-Flandrin distribution. 
%
%	X : signal (in time) to be analyzed. If X=[X1 X2], TFRDFLA 
%	   computes the cross-D-Flandrin distribution (Nx=length(X)).
%	T : time instant(s) on which the TFR is evaluated 
%	   					(default : 1:Nx).
%	FMIN,FMAX : respectively lower and upper frequency bounds of 
%	   the analyzed signal. These parameters fix the equivalent 
%	   frequency bandwidth (expressed in Hz). When unspecified, you
%	   have to enter them at the command line from the plot of the
%	   spectrum. FMIN and FMAX must be >0 and <=0.5.	 
%	N : number of analyzed voices (default : automatically determined).
%	TRACE : if nonzero, the progression of the algorithm is shown
%	                                 	(default : 0).
%	TFR : time-frequency matrix containing the coefficients of the 
%	   decomposition (abscissa correspond to uniformly sampled
%	   time, and ordonates correspond to a geometrically sampled 
%	   frequency). First row of TFR corresponds to the lowest 
%	   frequency. When called without output arguments, TFRDFLA
%	   runs TFRQVIEW.
%	F : vector of normalized frequencies (geometrically sampled 
%	   from FMIN to FMAX).
%
%	Example : 
%	 sig=altes(64,0.1,0.45); tfrdfla(sig); 
% 
%	See also all the time-frequency representations listed in
%	the file CONTENTS (TFR*).

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

global ratio_f

if (nargin == 0),
 error('At least one parameter required');
end;

[xrow,xcol] = size(X);
if nargin<=5, trace=0; end

if (nargin == 1),
 time=1:xrow; 
elseif (nargin==3),
 disp('FMIN will not be taken into account. Determine it with FMAX');
 disp('     from the following plot of the spectrum.'); 
elseif nargin==4,
 N=[];
end;

[trow,tcol] = size(time);
if (xcol==0)|(xcol>2),
 error('X must have one or two columns');
elseif (trow~=1),
 error('T must only have one row'); 
end; 

Mt = length(X); 

if trace, disp('D-Flandrin distribution'); end;

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

t = (1:Mt)-M-1;
T = xrow;

if nargin<=3,				        % fmin,fmax,N unspecified
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
 if isempty(fmin), fmin=fmindflt; end
 if isempty(fmax), fmax=fmaxdflt; end
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
ratio_f = fmax/fmin ;    
umax = fzero('umaxdfla',0); 
Teq = M/(fmax*umax);  
if Teq<2*M,
 M0 = round((2*M^2)/Teq-M)+1;
 T  = 2*(M+M0)-1;
else 
 M0 = 0;
end;
M1 = M+M0;

Nq= ceil((B*T*(1+2/R)*log((1+R/2)/(1-R/2)))/2);
Nmin = Nq-rem(Nq,2);
Ndflt = 2^nextpow2(Nmin);
if nargin<=3,
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

fmin_s = num2str(fmin) ; fmax_s = num2str(fmax) ; N_s = num2str(N) ;
if trace,
 disp(['frequency runs from ',fmin_s,' to ',fmax_s,' with ',N_s,' points']);
end


% Geometric sampling of the analyzed spectrum
k = 1:N;
q = (fmax/fmin)^(1/(N-1));
t = (1:Mt)-M-1;
geo_f  = fmin*(exp((k-1).*log(q)));
tfmatx = zeros(Mt,N);
tfmatx = exp(-2*i*t'*geo_f*pi);
S1 = s1'*tfmatx; 
S2 = s2'*tfmatx; 
clear tfmatx
S1(N+1:2*N) = zeros(1,N);
S2(N+1:2*N) = zeros(1,N);


% Mellin transform computation of the analyzed signal
p = 0:(2*N-1);
Mellin1 = fftshift(ifft(S1));
Mellin2 = fftshift(ifft(S2));
umin = -umax;
du = abs(umax-umin)/(2*M1);
u(1:2*M1) = umin:du:umax-du;
u(M1+1) = 0;
beta = (p/N-1)./(2*log(q));


% Computation of the Lambda(+/- u) dilations/compressions 
% of the analyzed signal
waf = zeros(2*M1,N);
for n = 1:2*M1,
 if trace, disprog(n,4*M1,10); end
 MX1 = exp(-(2*i*pi*beta+0.5)*2*log((1-u(n)/4))).*Mellin1;
 MX2 = exp(-(2*i*pi*beta+0.5)*2*log((1+u(n)/4))).*Mellin2;
 FX1 = fft(fftshift(MX1)) ;
 FX1 = FX1(1:N) ;
 FX2 = fft(fftshift(MX2)) ;
 FX2 = FX2(1:N) ;
 waf(n,:) = FX1.*conj(FX2);
end
waf = [waf(M1+1:2*M1,:) ; waf(1:M1,:)].*geo_f(ones(2*M1,1),:);
tffr = ifft(waf);  
tffr = real(rot90([tffr(M1+1:2*M1,:) ; tffr(1:M1,:)],-1));


% Conversion from [t.f,f] to [t,f] using a 1-D interpolation
tfr  = zeros(N,tcol);
Ts2  = (Mt-1)/2 ;
gamma = linspace(-geo_f(N)*Ts2,geo_f(N)*Ts2,2*M1) ;
alpha = (0.6*N-1)/0.4;
for n = 1:N,
 if trace, disprog(n+alpha,N+alpha,10); end
 ind = find(gamma>=-geo_f(n)*Ts2 & gamma<=geo_f(n)*Ts2);
 x = gamma(ind);
 y = tffr(n,ind);
 xi = (time-Ts2-1)*geo_f(n);
 v=interp1(x,y,xi,'spline');
 [l,r]=size(v);
 if (r==1),
  v=v';
 end;
 tfr(n,:)=v;
 clear v 
end 


t = time;
f = geo_f';

% Normalization
SP1 = fft(hilbert(s1)); 
SP2 = fft(hilbert(s2)); 
indmin = 1+round(fmin*(tcol-2));
indmax = 1+round(fmax*(tcol-2));
SP1ana=SP1(indmin:indmax);
SP2ana=SP2(indmin:indmax);
tfr=tfr*(SP1ana'*SP2ana)/integ2d(tfr,t,f)/N;

clear ratio_f

if (nargout==0),
 tfrqview(real(tfr),hilbert(real(X)),t,'tfrdfla',N,f);
end;

