function [tfr,rtfr,hat] = tfrrmsc(x,t,N,f0T,trace,K);
%TFRRMSC Reassigned Morlet Scalogram time-frequency distribution.
%	[TFR,RTFR,HAT] = TFRRMSC(X,T,N,F0T,TRACE) 
%	computes the Morlet scalogram and its reassigned version.
%                  
%	X     : analysed signal
%	T     : the time instant(s)           (default : 1:length(X))
%	N     : number of frequency bins      (default : length(X))
%	F0T   : time-bandwidth product of the mother wavelet 
%					      (default : 2.5)) 
%	TRACE : if nonzero, the progression of the algorithm is shown
%                                             (default : 0).
%	TFR,  : time-frequency representation and its reassigned
%	RTFR    version. When called without output arguments, 
%	        TFRRMSC runs TFRQVIEW.
%	HAT   : Complex matrix of the reassignment vectors.
%
%	Example :
%	 sig=fmlin(64,0.1,0.4); tfrrmsc(sig,1:64,64,2.1,1);
%
%	See also all the time-frequency representations listed in
%	 the file CONTENTS (TFR*)

%	F. Auger, January, April 1996.
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
if (nargin == 1),
 t=1:xrow; N=xrow; f0T=2.5; trace=0; K= 0.001;
elseif (nargin == 2), 
 N=xrow; f0T=2.5; trace=0; K= 0.001;
elseif (nargin == 3),
 f0T=2.5; trace=0; K=0.001;
elseif (nargin == 4),
 trace = 0; K=0.001;
elseif (nargin == 5),
 K= 0.001;
end;

if (N<0), error('N must be greater than zero'); end;
[trow,tcol] = size(t);

if (xcol~=1),
 error('X must have only one column');
elseif (trow~=1),
 error('T must only have one row');
elseif (f0T<=0),
 error('F0T must be positive');
elseif (length(f0T)~=1),
 error('F0T must be a scalar');
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

tfr= zeros(N,tcol); tf2= zeros(N,tcol);
if trace, disp('Morlet Scalogram'); end;

M=ceil(f0T*N*sqrt(2.0*log(1/K))); 
tau = 0:M+round(N/2); 
pi2 = 2.0*pi;
hstar = exp(-(tau/(N*f0T)).^2 /2.0) .* exp(-j*pi2*tau/N);
Thstar = tau.*hstar;

for m=1:round(N/2)-1
 if trace, disprog(m,N/2,10); end;
 factor=sqrt(m/(f0T*N));
 for icol=1:tcol,
   ti= t(icol);
   tauneg=1:min([ceil(M/m),ti-1]);
   taupos=0:min([ceil(M/m),xrow-ti]);
   % positive frequencies
   tfr(  1+m,icol)= hstar(1+taupos*m)*x(ti+taupos);
   tf2(  1+m,icol)=Thstar(1+taupos*m)*x(ti+taupos);
   if length(tauneg) > 0,
    tfr(1+m,icol)=tfr(1+m,icol) + conj( hstar(1+tauneg*m))*x(ti-tauneg);
    tf2(1+m,icol)=tf2(1+m,icol) - conj(Thstar(1+tauneg*m))*x(ti-tauneg);
   end;
   % negative frequencies
   tfr(N+1-m,icol)=conj( hstar(1+taupos*m))*x(ti+taupos);
   tf2(N+1-m,icol)=conj(Thstar(1+taupos*m))*x(ti+taupos);
   if length(tauneg) > 0,
    tfr(N+1-m,icol)=tfr(N+1-m,icol) +  hstar(1+tauneg*m)*x(ti-tauneg);
    tf2(N+1-m,icol)=tf2(N+1-m,icol) - Thstar(1+tauneg*m)*x(ti-tauneg);
   end;
 end;
 tfr(  1+m,:)=factor*tfr(  1+m,:); tf2(  1+m,:)=factor*tf2(  1+m,:)/m;
 tfr(N+1-m,:)=factor*tfr(N+1-m,:); tf2(N+1-m,:)=factor*tf2(N+1-m,:)/m;
end;

m=round(N/2); 
factor=sqrt(m/(f0T*N));
if trace, disprog(m,N/2,10); end
for icol=1:tcol,
 ti= t(icol);
 tauneg=1:min([ceil(M/m),ti-1]);
 taupos=0:min([ceil(M/m),xrow-ti]);
 tfr(  1+m,icol)= hstar(1+taupos*m)*x(ti+taupos);
 tf2(  1+m,icol)=Thstar(1+taupos*m)*x(ti+taupos);
 if length(tauneg) > 0,
  tfr(1+m,icol)=tfr(1+m,icol) + conj( hstar(1+tauneg*m))*x(ti-tauneg);
  tf2(1+m,icol)=tf2(1+m,icol) - conj(Thstar(1+tauneg*m))*x(ti-tauneg);
 end;
end;

tfr(1+m,:)=factor*tfr(1+m,:); 
tf2(1+m,:)=factor*tf2(1+m,:)/m;
tfr=tfr(:); tf2=tf2(:);

avoid_warn=find(tfr~=0.0);
tf2(avoid_warn)=tf2(avoid_warn)./tfr(avoid_warn); 
tfr=abs(tfr).^2;

if trace, disp('reassignment :'); end;
tfr=reshape(tfr,N,tcol);
tf2=reshape(tf2,N,tcol);

rtfr= zeros(N,tcol); 
Ex=mean(abs(x(min(t):max(t))).^2); Threshold=1.0e-6*Ex;
factor=2.0*pi*N*f0T*f0T;
for icol=1:tcol,
 if trace, disprog(icol,tcol,10); end;
 for jcol=1:N,
  if tfr(jcol,icol)>Threshold,
   icolhat= icol + round(real(tf2(jcol,icol)/Dt));
   icolhat=min(max(icolhat,1),tcol);
   m=rem(jcol+round(N/2)-2,N)-round(N/2)+1;
   jcolhat= jcol + round(imag(m*m*tf2(jcol,icol)/factor));
   jcolhat=rem(rem(jcolhat-1,N)+N,N)+1;
   rtfr(jcolhat,icolhat)= rtfr(jcolhat,icolhat)+tfr(jcol,icol);
   tf2(jcol,icol)= jcolhat + 1j * icolhat;
  else
   tf2(jcol,icol)=(1+j)*inf;
   rtfr(jcol,icol)=rtfr(jcol,icol)+tfr(jcol,icol);
  end;
 end;
end;

if (nargout==0),
 TFTBcontinue=1;
 while (TFTBcontinue==1),
  choice=menu ('Choose the representation:',...
               'stop',...
               'Morlet scalogram',...
               'reassigned Morlet scalogram');
  if (choice==1), TFTBcontinue=0;
  elseif (choice==2), 
   tfrqview(tfr,x,t,'tfrmsc',f0T);
  elseif (choice==3),
   tfrqview(rtfr,x,t,'tfrrmsc',f0T);
  end;
 end;
elseif (nargout>2),
 hat=tf2;
end;


