function [tfr,rtfr,hat] = tfrrgab(x,t,N,Nh,trace,K);
%TFRRGAB Reassigned Gabor spectrogram time-frequency distribution.
%	[TFR,RTFR,HAT] = TFRRGAB(X,T,N,NH,TRACE,K) 
%	computes the Gabor spectrogram and its reassigned version.
%	This particular window (a Gaussian window) allows a 20 % faster
%	algorithm than the TFRRSP function.
%                  
%	X     : analysed signal
%	T     : the time instant(s)           (default : 1:length(X))
%	N     : number of frequency bins      (default : length(X))
%	NH    : length of the gaussian window (default : N/4))
%	TRACE : if nonzero, the progression of the algorithm is shown
%                                             (default : 0).
%	K     : value at both extremities     (default 0.001)
%	TFR,  : time-frequency representation and its reassigned
%	RTFR    version. When called without output arguments, 
%	        TFRRGAB runs TFRQVIEW.
%	HAT   : Complex matrix of the reassignment vectors.
%
%	Example :
%	 sig=fmlin(128,0.1,0.4); tfrrgab(sig,1:128,128,19,1);
%
%	See also all the time-frequency representations listed in
%	 the file CONTENTS (TFR*)

%	F. Auger, May-July 1994, July 1995. 
%       Copyright (c) 1996 by CNRS(France). 
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
if (nargin <= 2),
 N=xrow;
end;

hlength=floor(N/4);
hlength=hlength+1-rem(hlength,2);

if (nargin == 1),
 t=1:xrow; 
end;

if (nargin <= 3),
 Nh=hlength; trace=0; K=0.001;
elseif (nargin == 4),
 trace = 0; K=0.001;
elseif (nargin == 5),
 K= 0.001;
end;

if (N<0),
 error('N must be greater than zero');
end;
[trow,tcol] = size(t);
if (xcol~=1),
 error('X must have only one column');
elseif (trow~=1),
 error('T must only have one row'); 
elseif (2^nextpow2(N)~=N & nargin==6),
 fprintf('For a faster computation, N should be a power of two\n');
end; 

if (rem(Nh,2)==0), 
 error('Nh must be odd'); 
elseif length(Nh)~=1,
 error('Nh must be a scalar');
end;

Nh2=Nh-2;
TFTBcontinue=1;
while TFTBcontinue,
 Nh2=Nh2+2;
 h=tftb_window(Nh2,'gauss',K^((Nh2-1)^2 /(Nh-1)^2)); 
 TFTBcontinue=(h(Nh2)*(Nh2-1)>2*K);
end;

K=K^((Nh2-1)^2 /(Nh-1)^2); Nh=Nh2; Lh=(Nh-1)/2; 
h=h; Th=h.*[-Lh:Lh]';

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
if trace, disp('Gabor spectrogram'); end;

for icol=1:tcol,
 if trace, disprog(icol,tcol,10); end;
 ti= t(icol); 
 tau=-min([round(N/2)-1,Lh,ti-1]):min([round(N/2)-1,Lh,xrow-ti]);
 indices= rem(N+tau,N)+1;
 norm_h=norm(h(Lh+1+tau));
 tfr(indices,icol)=x(ti+tau).*conj( h(Lh+1+tau)) /norm_h;
 tf2(indices,icol)=x(ti+tau).*conj(Th(Lh+1+tau)) /norm_h;
end ;
tfr=fft(tfr); tf2=fft(tf2);
tfr=tfr(:); tf2=tf2(:);  tf3=tf3(:);
avoid_warn=find(tfr~=0.0);
tf3(avoid_warn)=round(imag(2*log(K)*N*tf2(avoid_warn)./tfr(avoid_warn)/(2.0*pi*Lh^2)));
tf2(avoid_warn)=round(real(tf2(avoid_warn)./tfr(avoid_warn)/Dt));
tfr=abs(tfr).^2;
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
   icolhat= icol + tf2(jcol,icol);
   icolhat=min(max(icolhat,1),tcol);
   jcolhat= jcol - tf3(jcol,icol);
   %while (jcolhat<1),jcolhat=jcolhat+N; end;
   %while (jcolhat>N),jcolhat=jcolhat-N; end;
   jcolhat=rem(rem(jcolhat-1,N)+N,N)+1;
   rtfr(jcolhat,icolhat)=rtfr(jcolhat,icolhat) + tfr(jcol,icol) ;
   tf2(jcol,icol)=jcolhat + j * icolhat;
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
               'Gabor spectrogram',...
               'reassigned Gabor spectrogram');
  if (choice==1), TFTBcontinue=0;
  elseif (choice==2), 
   Q=round(tcol*N/xrow);
   tfrqview(tfr,x,t,'tfrgabor',tcol,Q,h);
  elseif (choice==3),
   tfrqview(rtfr,x,t,'tfrrgab',Nh);
  end;
 end;
elseif (nargout>2),
 hat=tf2;
end;

