function [waf,tau,theta] = ambifuwb(X,fmin,fmax,N,trace);
%AMBIFUWB Wide-band ambiguity function.
%       [WAF,TAU,THETA]=AMBIFUWB(X,FMIN,FMAX,N,TRACE) calculates
%	the asymetric wide-band ambiguity function.
%
%	X     : signal (in time) to be analyzed (the analytic associated
%	   signal is considered), of length Nx.
%       FMIN,FMAX : respectively lower and upper frequency bounds of
%	   the analyzed signal. When specified, these parameters fix
%	   the equivalent frequency bandwidth (in Hz). When unspecified, 
%	   you have to enter them at the command line from the plot of the
%	   spectrum. FMIN and FMAX must be >0 and <=0.5.
%       N     : number of Mellin points (default : automatically determined).
%	TRACE : if nonzero, the progression of the algorithm is shown
%					(default : 0).
%	AF    : matrix containing the coefficients of the ambiguity
%	   function. X-coordinate corresponds to the dual variable of 
%	   scale parameter ; Y-coordinate corresponds to time delay,
%	   dual variable of frequency.
%	   When called without output arguments, AMBIFUWB displays
%	   the squared modulus of the ambiguity function by means of
%	   contour.
%       TAU   : X-coordinate corresponding to time delay
%       THETA : Y-coordinate corresponding to the log(scale) variable
%
%
%	Example :
%        sig=altes(128,0.1,0.45); ambifuwb(sig);
%
%	See also AMBIFUNB.

%	P. Goncalvez - December 1995, O. Lemoine - August 1996.
%	Copyright (c) 1995 Rice University, CNRS (France)
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

[Nx,xcol] = size(X);
if (xcol==0)|(xcol>1),
 error('X must have one column');
end;
if nargin<=4, trace=0; end

if (nargin==2),
 disp('FMIN will not be taken into account. Determine it with FMAX');
 disp('     from the following plot of the spectrum.'); 
elseif (nargin==3),
 N=[];
end;

s = hilbert(real(X)); 
M = round(Nx/2);

t = (1:Nx)-M-1;
Tmin = 1;
Tmax = Nx;
T = Tmax-Tmin;

if nargin<=2,				        % fmin,fmax,N unspecified
 STF = fft(fftshift(s)); 
 sp  = (abs(STF(1:M))).^2; Maxsp=max(sp);
 f   = linspace(0,0.5,M+1) ; f=f(1:M);
 plot(f,sp) ; grid;
 xlabel('Normalized frequency');
 title('Analyzed signal energy spectrum');
 axis([0 1/2 0 1.2*Maxsp]) ; 
 indmin=min(find(sp>Maxsp/100));
 indmax=max(find(sp>Maxsp/100));
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
geo_f  = fmin*(exp((k-1).*log(q)));
tfmatx = zeros(Nx,N);
tfmatx = exp(-2*i*t'*geo_f*pi);
S = s.'*tfmatx; 
S = S(ones(1,Nx),:);
Sb = S.*tfmatx ;

tau = t;
S(:,N+1:2*N) = zeros(Nx,N);  S = S.';
Sb(:,N+1:2*N) = zeros(Nx,N); Sb = Sb.';

% Mellin transform computation of the analyzed signal
p=0:(2*N-1);
coef = exp(p/2.*log(q))'; 
MellinS = fftshift(ifft(S(:,1).*coef)).';
MellinS = MellinS(ones(1,Nx),:) ; MellinS = MellinS.';
for b=1:Nx,
 if trace, disprog(b,Nx,10); end
 MellinSb(:,b) = fftshift(ifft(Sb(:,b).*coef)) ;
end

k = 1:2*N;
beta = (p/N-1)/(2*log(q));

Scale = logspace(log10(fmin/fmax),log10(fmax/fmin),N) ;
waf = zeros(N,Nx) ;
MellinSSb = MellinS.*conj(MellinSb) ;

waf = ifft(MellinSSb,N);
No2=(N+rem(N,2))/2;
waf = [waf(No2+1:N,:) ; waf(1:No2,:)]; 

% Normalization
s=real(s);
SP = fft(hilbert(s)); 
indmin = 1+round(fmin*(Nx-2));
indmax = 1+round(fmax*(Nx-2));
SPana = SP(indmin:indmax);

waf=waf*norm(SPana)^2/waf(No2,M)/N;


theta = log(Scale) ;


if (nargout==0),
 contour(tau,theta,abs(waf.^2)); 
 grid
 if (exist('program_name')~=0),
  if (program_name == 'octave')
   xlabel('Delay'); ylabel('Doppler');
  end
 else
  xlabel('Delay'); ylabel('Doppler'); eval('shading interp');
 end
 title('Wide-band ambiguity function');
end;
