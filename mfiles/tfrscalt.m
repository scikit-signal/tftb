%function tfrscalt
%TFRSCALT Unit test for the time-frequency representation TFRSCALO.

%	O. Lemoine - June 1996. 

% We test each property of the corresponding TFR :

N=128;

% Covariance by translation in time 
t1=60; t2=70; f=0.3; W=0; 
sig1=amgauss(N,t1).*fmconst(N,f,t1); 
sig2=amgauss(N,t2).*fmconst(N,f,t2); 
tfr1=tfrscalo(sig1,1:N,W,0.1,0.4,128);  
tfr2=tfrscalo(sig2,1:N,W,0.1,0.4,128);        
[tr,tc]=size(tfr1);
nu=round(f*(tc-1)*2)+1;
tfr=tfr1-tfr2(:,modulo((1:tc)-t1+t2,tc));
if any(any(abs(tfr)>sqrt(eps))),
 error('tfrscalo test 1 failed');
end


% Covariance by dilation
t=N/2; f=0.2; T=2*sqrt(N); a=2; W=8; 
sig1=amgauss(N,t,T).*fmconst(N,f,t);
sig2=amgauss(a*N,a*t,T*a).*fmconst(a*N,f/a,a*t);
[tfr1,t1,f1]=tfrscalo(sig1,1:N  ,W,0.01,0.49,N);  
[tfr2,t2,f2]=tfrscalo(sig2,1:a*N,W,0.01,0.49,N);        
Max1=max(max(tfr1)); Max2=max(max(tfr2));
[I1,J1]=find(tfr1==Max1); [I2,J2]=find(tfr2==Max2);  
if abs(f1(I1)-a*f2(I2))>1e-2 | J2~=a*J1,
 error('tfrscalo test 2 failed');
end


% Reality of the TFR
sig=noisecg(N); W=5;
tfr=tfrscalo(sig,1:N,W,0.01,0.5,N);
if sum(any(abs(imag(tfr))>sqrt(eps)))~=0,
 error('tfrscalo test 3 failed');
end


% Energy conservation
sig=fmsin(N,.1,.4); W=6; Nf=2*N ;
[tfr,t,f]=tfrscalo(sig,1:N,W,0.01,0.49,Nf);
Es=norm(sig)^2/Nf;
Etfr=integ2d(tfr,t,f)/N;
if abs(Es-Etfr)>sqrt(eps),
 error('tfrscalo test 4 failed');
end


% Positivity
if any(any(tfr<0)),
 error('tfrscalo test 5 failed');
end


% Same energy in the time-scale plane for 2 gaussian atoms at different scales
sig=amgauss(256).*(fmconst(256,.15)+fmconst(256,.35));
[tfr,t,f]=tfrscalo(sig,1:256,12,.01,.49,512);
int1=integ2d(tfr(1:430,:),t,f(1:430));
int2=integ2d(tfr(431:512,:),t,f(431:512));
if abs(int1-int2)>1e-4,
 error('tfrscalo test 6 failed');
end


N=127;

% Covariance by dilation
t=ceil(N/2); f=0.2; T=2*sqrt(N); a=2; W=8; 
sig1=amgauss(N,t,T).*fmconst(N,f,t);
sig2=amgauss(a*N,a*t,T*a).*fmconst(a*N,f/a,a*t);
[tfr1,t1,f1]=tfrscalo(sig1,1:N  ,W,0.01,0.49,N);  
[tfr2,t2,f2]=tfrscalo(sig2,1:a*N,W,0.01,0.49,N);        
Max1=max(max(tfr1)); Max2=max(max(tfr2));
[I1,J1]=find(tfr1==Max1); [I2,J2]=find(tfr2==Max2);  
if abs(f1(I1(1))-a*f2(I2))>1e-2 | J2~=a*J1(1),
 error('tfrscalo test 7 failed');
end


% Reality of the TFR
sig=noisecg(N); W=5;
tfr=tfrscalo(sig,1:N,W,0.01,0.5,N);
if sum(any(abs(imag(tfr))>sqrt(eps)))~=0,
 error('tfrscalo test 8 failed');
end


% Energy conservation
sig=fmsin(N,.1,.4); W=6; Nf=2*N+1 ;
[tfr,t,f]=tfrscalo(sig,1:N,W,0.01,0.49,Nf);
SP = fft(hilbert(real(sig))); 
indmin = 1+round(0.01*(N-2));
indmax = 1+round(0.49*(N-2));
SPana = SP(indmin:indmax);
Es=SPana'*SPana/Nf;
Etfr=integ2d(tfr,t,f);
if abs(Es-Etfr)>1e-2,
 error('tfrscalo test 9 failed');
end


% Positivity
if any(any(tfr<0)),
 error('tfrscalo test 10 failed');
end


% Same energy in the time-scale plane for 2 gaussian atoms at different scales
sig=amgauss(N).*(fmconst(N,.15)+fmconst(N,.35));
[tfr,t,f]=tfrscalo(sig,1:N,12,.01,.49,2*N+1);
int1=integ2d(tfr(1:210,:),t,f(1:210));
int2=integ2d(tfr(211:255,:),t,f(211:255));
if abs(int1-int2)>5e-4,
 error('tfrscalo test 11 failed');
end

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
