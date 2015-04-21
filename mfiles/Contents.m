% Time-Frequency Toolbox.
% Version 1.0 January 1996
% Copyright (c) 1994-96 by CNRS (France) - RICE University (USA).
%
% Signal Generation Files
%
%   sigmerge - Add two signals with given energy ratio in dB.
%
%  Choice of the Instantaneous Amplitude 
%   amexpo1s - Generate one-sided exponential amplitude modulation.
%   amexpo2s - Generate bilateral exponential amplitude modulation.
%   amgauss  - Generate gaussian amplitude modulation.
%   amrect   - Generate rectangular amplitude modulation.
%   amtriang - Generate triangular amplitude modulation.
%   
%  Choice of the Instantaneous Frequency
%   fmconst  - Signal with constant frequency modulation.
%   fmhyp    - Signal with hyperbolic frequency modulation.
%   fmlin    - Signal with linear frequency modulation.
%   fmodany  - Signal with arbitrary frequency modulation.
%   fmpar    - Signal with parabolic frequency modulation.
%   fmpower  - Signal with power-law frequency modulation.
%   fmsin    - Signal with sinusoidal frequency modulation.
%   gdpower  - Signal with a power-law group delay.
%
%  Choice of Particular Signals
%   altes    - Altes signal in time domain.
%   anaask   - Amplitude Shift Keyed (ASK) signal.
%   anabpsk  - Binary Phase Shift Keyed (BPSK) signal.
%   anafsk   - Frequency Shift Keyed (FSK) signal.
%   anapulse - Analytic projection of unit amplitude impulse signal.
%   anaqpsk  - Quaternary Phase Shift Keyed (QPSK) signal.
%   anasing  - Lipschitz singularity.
%   anastep  - Analytic projection of unit step signal.
%   atoms    - Linear combination of elementary Gaussian wave packets. 
%   dopnoise - Generate complex Doppler random signal.
%   doppler  - Generate complex Doppler signal.
%   klauder  - Klauder wavelet in time domain.
%   mexhat   - Mexican hat wavelet in time domain.
%   tftb_window - Window generation (previously window.m).
%
%  Addition of Noise
%   noisecg  - Analytic complex gaussian noise.
%   noisecu  - Analytic complex uniform noise.
%
%  Modification
%   scale    - Scale signal using Mellin transform.
%   
%
% Processing Files
%
%  Time-Domain Processing
%   ifestar2 - Instantaneous frequency estimation using AR2 modelisation.
%   instfreq - Instantaneous frequency estimation.
%   loctime  - Time localization characteristics.
%	
%  Frequency-Domain Processing
%   fmt      - Fast Mellin transform
%   ifmt     - Inverse fast Mellin transform.
%   locfreq  - Frequency localization characteristics.
%   sgrpdlay - Group delay estimation.
%
%  Linear Time-Frequency Processing
%   tfrgabor - Gabor representation.
%   tfrstft  - Short time Fourier transform.
%
%  Bilinear Time-Frequency Processing in the Cohen's Class
%   tfrbj    - Born-Jordan distribution.                     
%   tfrbud   - Butterworth distribution.                   
%   tfrcw    - Choi-Williams distribution.                   
%   tfrgrd   - Generalized rectangular distribution.        
%   tfrmh    - Margenau-Hill distribution.                    
%   tfrmhs   - Margenau-Hill-Spectrogram distribution.     
%   tfrmmce  - MMCE combination of spectrograms. 
%   tfrpage  - Page distribution.                          
%   tfrpmh   - Pseudo Margenau-Hill distribution.          
%   tfrppage - Pseudo Page distribution.                 
%   tfrpwv   - Pseudo Wigner-Ville distribution.            
%   tfrri    - Rihaczek distribution.                       
%   tfrridb  - Reduced interference distribution (Bessel window).
%   tfrridh  - Reduced interference distribution (Hanning window).
%   tfrridn  - Reduced interference distribution (binomial window).
%   tfrridt  - Reduced interference distribution (triangular window).
%   tfrsp    - Spectrogram.                     
%   tfrspwv  - Smoothed Pseudo Wigner-Ville distribution.   
%   tfrwv    - Wigner-Ville distribution.
%   tfrzam   - Zhao-Atlas-Marks distribution.               
%
%  Bilinear Time-Frequency Processing in the Affine Class
%   tfrbert  - Unitary Bertrand distribution.
%   tfrdfla  - D-Flandrin distribution.
%   tfrscalo - Scalogram, for Morlet or Mexican hat wavelet.
%   tfrspaw  - Smoothed Pseudo Affine Wigner distributions.
%   tfrunter - Unterberger distribution, active or passive form.
%   
%  Reassigned Time-Frequency Processing
%   tfrrgab  - Reassigned Gabor spectrogram.
%   tfrrmsc  - Reassigned Morlet Scalogram time-frequency distribution.
%   tfrrpmh  - Reassigned Pseudo Margenau-Hill distribution.
%   tfrrppag - Reassigned Pseudo Page distribution.
%   tfrrpwv  - Reassigned Pseudo Wigner-Ville distribution.
%   tfrrsp   - Reassigned Spectrogram. 
%   tfrrspwv - Reassigned Smoothed Pseudo WV distribution.
%
%  Ambiguity Functions
%   ambifunb - Narrow band ambiguity function.
%   ambifuwb - Wide band ambiguity function.
%
%  Post-Processing or Help to the Interpretation
%   friedman - Instantaneous frequency density.
%   holder   - Estimate the Holder exponent through an affine TFR.
%   htl      - Hough transform for detection of lines in images.
%   margtfr  - Marginals and energy of a time-frequency representation.
%   midscomp - Mid-point construction used in the interference diagram. 
%   momftfr  - Frequency moments of a time-frequency representation.
%   momttfr  - Time moments of a time-frequency representation.
%   plotsid  - Schematic interference diagram of FM signals. 
%   renyi    - Measure Renyi information.
%   ridges   - Extraction of ridges.
%   tfrideal - Ideal TFR for given frequency laws.
%
%  Visualization & backup
%   plotifl  - Plot normalized instantaneous frequency laws.
%   tfrqview - Quick visualization of time-frequency representations.
%   tfrview  - Visualization of time-frequency representations.
%   tfrparam - Return the paramaters needed to display (or save) a TFR.
%   tfrsave  - Save the parameters of a time-frequency representation.
%
%
% Other 
%   disprog  - Display progression of a loop.
%   divider  - Find dividers of integer such that product equals integer. 
%   dwindow  - Derive a window.
%   integ    - Approximate integral.
%   integ2d  - Approximate 2-D integral.
%   izak     - Inverse Zak transform.
%   istfr1   - returns true if a distribution is Cohen's class type 1 (spectrogram)
%   istfr2   - returns true if a distribution is Cohen's class type 2 (Wigner-Ville)
%   istfraff - returns true is a distribution is an affine class member
%   kaytth   - Kay-Tretter filter computation.
%   modulo   - Congruence of a vector.
%   movcw4at - Four atoms rotating, analyzed by the Choi-Williams distribution.
%   movpwdph - Influence of a phase-shift on the interferences of the pWVD.
%   movpwjph - Influence of a jump of phase on the interferences of the pWVD.
%   movsc2wv - Movie illustrating the passage from the scalogram to the WVD.
%   movsp2wv - Movie illustrating the passage from the spectrogram to the WVD.
%   movwv2at - Oscillating structure of the interferences of the WVD.  
%   odd      - Round towards nearest odd value.
%   sigmerge - Add two signals with given energy ratio in dB.
%   zak	     - Zak transform.

%   griffitc - Test signal example C of Griffiths' paper.
%   lambdak  - Evaluate lambda function for Affine Wigner distribution.
%   umaxbert - Determination of the maximum value of u for Bertrand distribution.
%   umaxdfla - Determination of the maximum value of u for D-Flandrin distribution.
%   umaxunte - Determination of the maximum value of u for Unterberger distribution.

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
