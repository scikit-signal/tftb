function g = gaussn(f0,n)

% function gn = gaussn(f0,n) : generates the order n derivative of the 
% gaussian window, centered at frequency f0
% The wavelet gn is real, but it is its analytic form that is 
% synthetized first. The real part is then normalized by the analytical
% value of its energy (so, L2 normalization!)

% alpha : std of the initial Gaussian 
% g_0(t) = sqrt(2) (alpha/2*pi)^(1/4) exp(-alpha t^2)

alpha = 2*pi^2*f0^2 / n ;

% N = number of points, depends on "n". The support of the wavelet 
% at tolerance 10^(-3) is determined experimentaly and modeled by
% a 4th order polynome (with reference frequency f0 is 0.1) 

PolyCoef = [-0.00000001420054   0.00001199042535  -0.00403466080616 ...
      0.92639865727446 20.64720217778273] ;
N = ceil((0.1/f0).*polyval(PolyCoef,n)) ;

f = linspace(0,1,2*N+1) ; f = f(1:2*N) ;

% Synthesis of the "analytic" wavelet Gn in the frequency domain 

G = (2*pi/alpha)^(1/4) * (2*i*pi*f).^n .* exp(-(pi^2*f.^2)./alpha) ;

% Calculus of the wavelet energy (theoretical expression, cf. Gratshteyn,
% sec. 3.461, form. 2)

p = (2*pi^2)/alpha ;
E = (2*pi/alpha)^(1/2) * (2*pi).^(2*n) * ...
    (fact2(2*n-1))/(2*(2*p)^n)*sqrt(pi/p) ;
if isinf(E) | isnan(E) | E == 0
  E = integ(abs(G).^2,f) ;
end

% L2 Normalized wavelet in time domain 

g = real ( fftshift(ifft(G))./(sqrt(E/2)) ) ;
g = g(2:2*N) ;
t = -(N-1):(N-1) ;

% subplot(211) ;
% plot(f,abs(G)) ; title ('Wavelet Spectrum') ;
% xlabel('Frequency') ; grid 
% subplot(212)
% plot(t,g)
% title('Wavelet in time')
% xlabel('Time') ; grid


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
