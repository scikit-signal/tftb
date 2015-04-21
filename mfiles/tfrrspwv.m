function [tfr,rtfr,hat] = tfrrspwv(x,t,N,g,h,trace);
%TFRSPWV Reassigned smoothed pseudo Wigner-Ville distribution.
%	[TFR,RTFR,HAT] = TFRRSPWV(X,T,N,G,H,TRACE) 
%	computes the smoothed pseudo Wigner-Ville distribution and its
%	reassigned version.
% 
%	X     : analysed signal.
%	T     : the time instant(s)      (default : 1:length(X)).
%	N     : number of frequency bins (default : length(X)).
%	G     : time smoothing window, G(0) being forced to 1. 
%	                                 (default : Hamming(N/10)). 
%	H     : frequency smoothing window, H(0) being forced to 1
%	                                 (default : Hamming(N/4)).
%	TRACE : if nonzero, the progression of the algorithm is shown
%	                                 (default : 0).
%	TFR,  : time-frequency representation and its reassigned
%	RTFR  : version. When called without output arguments, 
%	         TFRRSPWV runs TFRQVIEW.
%	HAT   : Complex matrix of the reassignment vectors.
%
%	Example :
%	 sig=fmlin(128,0.05,0.15)+fmlin(128,0.3,0.4); t=1:2:128; 
%	 g=tftb_window(15,'Kaiser'); h=tftb_window(63,'Kaiser');  
%	 tfrrspwv(sig,t,64,g,h,1);
%
%	See also all the time-frequency representations listed in
%	 the file CONTENTS (TFR*)

%	F. Auger, May-July 1994, July 1995.
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

if (nargin == 0),
 error('At least 1 parameter required');
end;
[xrow,xcol] = size(x);
if (xcol~=1),
 error('X must have only one column');
end;

if (nargin <= 2),
 N=xrow;
elseif (N<0),
 error('N must be greater than zero');
elseif (2^nextpow2(N)~=N),
 fprintf('For a faster computation, N should be a power of two\n');
end;

hlength=floor(N/4); hlength=hlength+1-rem(hlength,2); 
glength=floor(N/10);glength=glength+1-rem(glength,2);

if (nargin == 1),
 t=1:xrow; g = tftb_window(glength); h = tftb_window(hlength); trace = 0;
elseif (nargin == 2)|(nargin == 3),
 g = tftb_window(glength); h = tftb_window(hlength); trace = 0;
elseif (nargin == 4),
 h = tftb_window(hlength); trace = 0;
elseif (nargin == 5),
 trace = 0;
end;

[trow,tcol] = size(t);
if (trow~=1),
 error('T must only have one row'); 
end; 

[grow,gcol]=size(g); Lg=(grow-1)/2; % g=g/sum(g);
if (gcol~=1)|(rem(grow,2)==0),
 error('G must be a smoothing window with odd length'); 
end;

[hrow,hcol]=size(h); Lh=(hrow-1)/2; h=h/h(Lh+1);
if (hcol~=1)|(rem(hrow,2)==0),
  error('H must be a smoothing window with odd length');
end;

if (tcol==1),
 Dt=1; 
else
 Deltat=t(2:tcol)-t(1:tcol-1); 
 Mini=min(Deltat); Maxi=max(Deltat);
 if (Mini~=Maxi),
  error('The time instants must be regularly sampled.');
 else
  Dt=Mini;
 end;
 clear Deltat Mini Maxi;
end;

tfr= zeros(N,tcol); tf2= zeros(N,tcol); tf3= zeros(N,tcol);
if trace, disp('Smoothed pseudo Wigner-Ville distribution'); end;
Dh=dwindow(h); % Tg=g.*[-Lg:Lg]'; 
for icol=1:tcol,
 ti= t(icol);
 taumax=min([ti+Lg-1,xrow-ti+Lg,round(N/2)-1,Lh]);
 if trace, disprog(icol,tcol,10); end;

 % tau=0
 points= -min([Lg,xrow-ti]):min([Lg,ti-1]);
 g2=g(Lg+1+points); g2=g2/sum(g2); Tg2= g2 .* points.' ;
 xx= x(ti-points) .* conj(x(ti-points));
 tfr(1,icol)= sum( g2 .* xx) ; 
 tf2(1,icol)= sum( Tg2 .* xx) ;
 tf3(1,icol)= Dh(Lh+1) * tfr(1,icol) ;

 for tau=1:taumax,
  points= -min([Lg,xrow-ti-tau]):min([Lg,ti-tau-1]);
  g2=g(Lg+1+points); g2=g2/sum(g2); Tg2= g2 .* points.' ;
  xx=x(ti+tau-points,1) .* conj(x(ti-tau-points));
  tfr(  1+tau,icol)= sum( g2 .* xx);
  tf3(  1+tau,icol)=Dh(Lh+tau+1) * tfr(  1+tau,icol) ;
  tfr(  1+tau,icol)= h(Lh+tau+1) * tfr(  1+tau,icol) ;
  tf2(  1+tau,icol)= h(Lh+tau+1) * sum(Tg2 .* xx);

  tfr(N+1-tau,icol)= sum( g2 .* conj(xx));
  tf3(N+1-tau,icol)=Dh(Lh-tau+1) * tfr(N+1-tau,icol);
  tfr(N+1-tau,icol)= h(Lh-tau+1) * tfr(N+1-tau,icol);
  tf2(N+1-tau,icol)= h(Lh-tau+1) * sum(Tg2 .* conj(xx));
 end;
end;

tfr=real(fft(tfr));
tf2=real(fft(tf2));
tf3=imag(fft(tf3));
tfr=tfr(:); tf2=tf2(:); tf3=tf3(:);
avoid_warn=find(tfr~=0);
tf2(avoid_warn)=round(tf2(avoid_warn)./tfr(avoid_warn)/Dt);
tf3(avoid_warn)=round(N*tf3(avoid_warn)./tfr(avoid_warn)/(2.0*pi));
if trace, fprintf ('\nreassignment: \n'); end;
tfr=reshape(tfr,N,tcol);
tf2=reshape(tf2,N,tcol);
tf3=reshape(tf3,N,tcol);


rtfr= zeros(N,tcol); 
Ex=mean(abs(x(min(t):max(t))).^2); Threshold=1.0e-6*Ex;
for icol=1:tcol,
 if trace, disprog(icol,tcol,10); end;
 for jcol=1:N,
  if abs(tfr(jcol,icol))>Threshold,
   icolhat= icol - tf2(jcol,icol);
   icolhat=min(max(icolhat,1),tcol);
   jcolhat= jcol - tf3(jcol,icol);
   jcolhat=rem(rem(jcolhat-1,N)+N,N)+1;
   rtfr(jcolhat,icolhat)= rtfr(jcolhat,icolhat) + tfr(jcol,icol);
   tf2(jcol,icol)= jcolhat + j * icolhat;
  else 
   tf2(jcol,icol)=inf*(1+j);
   rtfr(jcol,icol)=rtfr(jcol,icol) + tfr(jcol,icol) ;
  end;
 end;
end;

if trace, fprintf('\n'); end;
clear tf3; 

if (nargout==0),
 TFTBcontinue=1;
 while (TFTBcontinue==1),
  choice=menu ('Choose the representation:',...
               'stop',...
               'smoothed pseudo Wigner-Ville distribution',...
               'reassigned smoothed pseudo Wigner-Ville distribution');
  if (choice==1), TFTBcontinue=0;
  elseif (choice==2), 
   tfrqview(tfr,x,t,'tfrspwv',g,h);
  elseif (choice==3),
   tfrqview(rtfr,x,t,'tfrrspwv',g,h);
  end;
 end;
elseif (nargout>2),
 hat=tf2;
end;
