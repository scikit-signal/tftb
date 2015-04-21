function [stan,rtfr,hat] = tfrrstan(x,t,N,G,h,trace);
%TFRRSTAN Reassigned Stankovic distribution.
%	[TFR,RTFR,HAT] = TFRRSTAN(X,T,N,H,TRACE) 
%	computes the Stankovic distribution and its reassigned version.
% 
%	X     : analysed signal.
%	T     : the time instant(s)      (default : 1:length(X)).
%	N     : number of frequency bins (default : length(X)).
%       G     : frequency averaging window
%	h     : stft window              (default : Hamming(N/4)).
%	TRACE : if nonzero, the progression of the algorithm is shown
%	                                 (default : 0).
%	TFR,  : time-frequency representation and its reassigned
%	RTFR    version. When called without output arguments, 
%	        TFRRSTAN runs TFRQVIEW.
%	HAT   : Complex matrix of the reassignment vectors.
%
%	Example :
%	 sig=fmlin(128,0.1,0.4); t=1:2:128;
%	 G=tftb_window(9,'hanning'); h=tftb_window(61,'hanning'); tfrrstan(sig,t,128,G,h,1);
%
%	See also  all the time-frequency representations listed in
%	 the file CONTENTS (TFR*)

%	F. Auger, August, September 1997.
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
 t=1:xrow; G=[0.25; 0.5; 0.25]; h = tftb_window(hlength); trace=0;
elseif (nargin == 2)|(nargin == 3),
 G=[0.25; 0.5; 0.25]; h = tftb_window(hlength); trace=0;
elseif (nargin == 4),
 h = tftb_window(hlength); trace = 0;
elseif (nargin == 5),
 trace = 0;
end;

[Grow,Gcol]=size(G); LG=(Grow-1)/2; 
if (Gcol~=1)|(rem(Grow,2)==0),
 error('G must be a smoothing window with odd length');
end;

if (N<0),
 error('N must be greater than zero');
end;
[trow,tcol] = size(t);
if (xcol~=1),
 error('X must have only one column');
elseif (trow~=1),
 error('T must only have one row'); 
elseif (2^nextpow2(N)~=N),
 fprintf('For a faster computation, N should be a power of two\n');
end; 

[hrow,hcol]=size(h); Lh=(hrow-1)/2; 
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

stan = zeros(N,tcol); 
rtfr = zeros(N,tcol); 
hat  = zeros(N,tcol);

if trace, disp('Stankovic distribution (with reassignement)'); end;
Ex=mean(abs(x(min(t):max(t))).^2); Threshold=1.0e-3*Ex;
Dh=dwindow(h); Th=h.*[-Lh:Lh]';
for icol=1:tcol,
 if trace, disprog(icol,tcol,10); end;
 ti= t(icol); 
 tau=-min([round(N/2)-1,Lh,ti-1]):min([round(N/2)-1,Lh,xrow-ti]);
 indices= rem(N+tau,N)+1;
 tfr= zeros(N,3); 
 tfr(indices,1)=x(ti+tau).*conj( h(Lh+1-tau));
 tfr(indices,2)=x(ti+tau).*conj(Th(Lh+1-tau));
 tfr(indices,3)=x(ti+tau).*conj(Dh(Lh+1-tau));
 tfr=fft(tfr); 
 stan(:,icol)=G(LG+1) * abs(tfr(:,1)).^2;

 for jcol=1:N,
  stanTh=G(LG+1)*tfr(jcol,1)*conj(tfr(jcol,2));
  stanDh=G(LG+1)*tfr(jcol,1)*conj(tfr(jcol,3));
  for kstan=1:min(N/2-1,LG),
   stanbefore=rem(rem(jcol-kstan-1,N)+N,N)+1;
   stanafter =rem(rem(jcol+kstan-1,N)+N,N)+1;
   stan(jcol,icol)= stan(jcol,icol) ...
                  + G(LG+1-kstan)*tfr(stanbefore,1)*conj(tfr(stanafter ,1)) ...
                  + G(LG+1+kstan)*tfr(stanafter ,1)*conj(tfr(stanbefore,1));
   stanTh= stanTh + G(LG+1-kstan)*tfr(stanbefore,1)*conj(tfr(stanafter ,2)) ...
                  + G(LG+1+kstan)*tfr(stanafter ,1)*conj(tfr(stanbefore,2));
   stanDh= stanDh + G(LG+1-kstan)*tfr(stanbefore,1)*conj(tfr(stanafter ,3)) ...
                  + G(LG+1+kstan)*tfr(stanafter ,1)*conj(tfr(stanbefore,3));
  end;
  stan(jcol,icol)=real(stan(jcol,icol));
  if abs(stan(jcol,icol))>Threshold,
   icolhat = round(icol - real(stanTh/stan(jcol,icol))/Dt);
   jcolhat = round(jcol - N*imag(stanDh/stan(jcol,icol)/(2.0*pi)));
   
   jcolhat= rem(rem(jcolhat-1,N)+N,N)+1;
   icolhat= min(max(icolhat,1),tcol);
   
   rtfr(jcolhat,icolhat)=rtfr(jcolhat,icolhat) + stan(jcol,icol) ;
   hat(jcol,icol)= jcolhat + j * icolhat;
   %fprintf('%12.3f %12.3f , %12.3f %12.3f \n',jcol,icol,jcolhat,icolhat);
  else
   rtfr(jcol,icol)=rtfr(jcol,icol) + stan(jcol,icol);
   hat(jcol,icol)= jcol + j * icol;
  end;
 end;

end ;

if trace, fprintf('\n'); end;

if (nargout==0),
 TFTBcontinue=1;
 while (TFTBcontinue==1),
  choice=menu ('Choose the representation:',...
               'stop',...
               'Stankovic distribution',...
               'reassigned Stankovic distribution');
  if (choice==1), TFTBcontinue=0;
  elseif (choice==2), 
   tfrqview(stan,x,t,'type1');
  elseif (choice==3),
   tfrqview(rtfr,x,t,'type1');
  end;
 end;
end;
