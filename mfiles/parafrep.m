function [spec,freqs,A]=parafrep(Rx,N,method,q);
% [spec,freqs]=parafrep(Rx,N,method,q) parametric frequency representation
% of a signal.
%
% Rx     : correlation matrix of size (p+1)x(p+1)
% N      : number of frequency bins between 0 and 0.5
% method : can be either 'ar', 'periodogram', 'capon', 'capnorm', 'lagunas',
%          or 'genlag'.
% q      : parameter for the generalized Lagunas method.
%
% noise=rand(1000,1); signal=filter([1 0 0],[1 1 1],noise);
% figure(1);parafrep(correlmx(signal,2,'hermitian'),128,'AR');title('AR (2)');
% figure(2);parafrep(correlmx(signal,4,'hermitian'),128,'Capon');title('Capon (4)');
% figure(3);parafrep(correlmx(signal,2,'hermitian'),128,'lagunas');title('Lagunas (2)');
% figure(4);parafrep(correlmx(signal,40,'hermitian'),128,'periodogram');title('periodogram (40)');

% F. Auger, july 1998, april 99.

if nargin<1,
 error('At least one parameter required');
elseif nargin==1,
 N=128; method='capon';
elseif nargin==2,
 method='capon';
end;

[Rxrow,Rxcol]=size(Rx);
if (Rxrow ~= Rxcol),
 error('Rx must be a square matrix');
end;

p=Rxrow-1;
freqs=linspace(0,0.5,N);
spec=zeros(N,1);

method=upper(method);
if strcmp(method,'AR'),
 Un=ones(Rxrow,1); Rxm1Un= (Rx\Un); 
 P1=real(Un'*Rxm1Un); A=Rxm1Un/P1; 
 for ifreq=1:N, 
  Z=exp(2j*pi*freqs(ifreq)*(0:Rxrow-1)'); 
  spec(ifreq)=P1 ./ abs(Z' * A)^2 ;
 end;
elseif strcmp(method,'PERIODOGRAM'),
 for ifreq=1:N, 
  Z=exp(2j*pi*freqs(ifreq)*(0:Rxrow-1)'); 
  spec(ifreq)=real(Z' * Rx *Z)/(p+1)^2;
 end; 
elseif strcmp(method,'CAPON'),
 for ifreq=1:N, 
  Z=exp(2j*pi*freqs(ifreq)*(0:Rxrow-1)'); 
  spec(ifreq)=1.0 / real(Z' * (Rx\Z));
 end; 
elseif strcmp(method,'CAPNORM'),
 for ifreq=1:N, 
  Z=exp(2j*pi*freqs(ifreq)*(0:Rxrow-1)'); 
  spec(ifreq)=(p+1) / real(Z' * (Rx\Z));
 end; 
elseif strcmp(method,'LAGUNAS'),
 for ifreq=1:N, 
  Z=exp(2j*pi*freqs(ifreq)*(0:Rxrow-1)'); 
  Rxm1Z=Rx\Z; spec(ifreq)=real(Z' * Rxm1Z)/real(Z' * (Rx\Rxm1Z));
 end; 
elseif strcmp(method,'GENLAG'),
 for ifreq=1:N, 
  Z=exp(2j*pi*freqs(ifreq)*(0:Rxrow-1)'); 
  Rxqm1Z=(Rx)^q \Z; spec(ifreq)=real(Z' * Rx * Rxqm1Z)/real(Z' * (Rx\Rxqm1Z));
 end; 
else
 error('unknown frequency representation');
end;

if (nargout==0),
 figure(gcf); plot(freqs,10.0*log10(spec)); grid;
 xlabel('normalized frequency');
 ylabel('DSP  (dB)');
end;


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
