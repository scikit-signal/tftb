=================================================
Gallery of Time-Frequency Distributions in PyTFTB
=================================================

.. toctree::
  :maxdepth: 2

This is a gallery of various signals that are generated, processed, and
visualized through PyTFTB. Generation of signals falls broadly under three
categories:

- Frequency modulated signals
- Amplitude modulated signals
- Analytic signals
- Miscellaneous generation utilities (noise generation etc.)

The following plots display some of these signals and their different time-frequency
representations.

Linear frequency modulation (chirp)
-----------------------------------
  
.. plot:: _gallery/plot_1_3_1_chirp.py

The energy spectrum of the chirp.

.. plot:: _gallery/plot_1_3_1_chirp_spectrum.py

The Wigner-Ville distribution of the chirp.

.. plot:: _gallery/plot_1_3_1_chirp_wv.py

Noisy chirp
-----------

A chirp with white Gaussian noise.

.. plot:: _gallery/plot_1_3_1_noisy_chirp.py

Enery spectrum of the noisy chirp.

.. plot:: _gallery/plot_1_3_1_noisy_chirp_spectrum.py

The Wigner-Ville distribution of the noisy chirp.

.. plot:: _gallery/plot_1_3_1_noisy_chirp_wv.py

Localization of time-frequency characteristics
----------------------------------------------

A Gaussian amplitude modulated sinusoid with time center at t = 127 and a
frequency center at f = 0.35

.. plot:: _gallery/plot_2_2_1_time_freq_localization.py

Instantaneous frequency and group delay estimation
--------------------------------------------------

The estimated instantaneous frequency of a chirp.

.. plot:: _gallery/plot_2_3_instantaneous_frequency.py

The estimated group delay of a chirp.

.. plot:: _gallery/plot_2_4_group_delay.py

Comparison of the group delay and instantaneous frequencies of two signals
having high and low time-bandwidth product respectively.

.. plot:: _gallery/plot_2_4_grp_delay_inst_freq_comprison.py

Synthesis of monocomponent, nonstationary signals
---------------------------------------------------

A monocomponent, nonstationary signal with linear frequency modulation and
Gaussian amplitude modulation.

.. plot:: _gallery/plot_2_6_monocomp_nonstat_linfreq_gaussamp.py

A monocomponent, nonstationary signal with constant frequency modulation and
exponential amplitude modulation.

.. plot:: _gallery/plot_2_6_monocomp_nonstat_constfreq_expamp.py

A doppler wave.

.. plot:: _gallery/plot_2_6_monocomp_nonstat_doppler.py

Gaussian transient signal in -10 dB colored Gaussian noise.

.. plot:: _gallery/plot_2_6_monocomp_nonstat_colored_gaussian_noise.py

Synthesis of multicomponent, nonstationary signals
--------------------------------------------------

The instantaneous frenquency and group delay estimation of a multi-component
signal.

.. plot:: _gallery/plot_2_7_multicomp_nonstat_instfreq_grpdlay.py

Short-time Fourier transform of the multicomponent nonstationary signal.

.. plot:: _gallery/plot_2_7_multicomp_nonstat_stft.py

Examples of short-time Fourier transform
----------------------------------------

An audio signal and its energy spectrum:

.. plot:: _gallery/plot_3_1_2_spectrum.py

Short time Fourier transform of the audio signal.

.. plot:: _gallery/plot_3_1_2_stft.py

Patterns corresponding to both the sound in the signal, and to the harmonics
can be see.

Perfect time resolution in STFT:

.. plot:: _gallery/plot_3_1_4_time_resolution.py

Perfect frequency resolution in STFT:

.. plot:: _gallery/plot_3_1_4_frequency_resolution.py

STFT with a long smoothing window, showing good resolution in frequency, but
poor resolution in time:

.. plot:: _gallery/plot_3_1_4_atoms_hamming_stft.py

STFT with a short smoothing window, showing good resoultion in time, but
slightly worse resolution in freqeuncy:

.. plot:: _gallery/plot_3_1_4_atoms_short_hamming_stft.py

The Gabor Representation
------------------------

The biorthonormal window used by the Gabor representation, at the critical
sampling rate.

.. plot:: _gallery/plot_3_3_2_biorthonormal_window.py

The Gabor representation of a chirp, at the critical sampling rate.

.. plot:: _gallery/plot_3_3_2_biorthonormal_window_gabor.py

The Gabor representation of the same chirp, but highly oversampled:

.. plot:: _gallery/plot_3_3_2_biorthonormal_window_oversampled.py

The Spectrogram Representation
------------------------------

Spectrogram of two parallel chirps, with a short Gaussian smoothing window.
Some ambiguity in frequency can be seen.

.. plot:: _gallery/plot_3_4_1_chirps_spectrogram_short_gaussian.py

Spectrogram of two parallel chirps, with a long Gaussian smoothing window. Here
the ambiguity is in time.

.. plot:: _gallery/plot_3_4_1_chirps_spectrogram_long_gaussian.py

Spectrogram of two chirps, more distant than in the previous example, with a
short Gaussian smoothing window.

.. plot:: _gallery/plot_3_4_1_distant_components_short_gaussian.py

Spectrogram of the two distant chirps, with a long Gaussian smoothing window.

.. plot:: _gallery/plot_3_4_1_distant_components_long_gaussian.py

The Scalogram Representation
----------------------------

The Morlet scalogram of a Dirac impulse.

.. plot:: _gallery/plot_3_4_2_morlet_scalogram_dirac_impulse.py

The Morlet scalogram of two complex simulatneous sinusoids.

.. plot:: _gallery/plot_3_4_2_morlet_scalogram_complex_sinusoids.py

The Wigner-Ville Distribution
-----------------------------

The Wigner-Ville distribution of a linear chirp:

.. plot:: _gallery/plot_4_1_1_wv_wireframe.py

The Cohen's class of time-frequency representations
---------------------------------------------------

WVD of a signal containing a Gaussian atom and a complex sinusoid.

.. plot:: _gallery/plot_4_1_2_sin_gauss_wv.py

Pseudo WV distribution of the previous signal.

.. plot:: _gallery/plot_4_1_2_sin_gauss_pwv.py

WVD of chirps with Gaussian amplitudes and different slopes.

.. plot:: _gallery/plot_4_1_3_chirps_wvd.py

Narrow band ambiguity function of the previous signal.

.. plot:: _gallery/plot_4_1_3_chirps_ambifunb.py

The Margenau-Hill representation
--------------------------------

Margenau-Hill representation of two Gaussian atoms.

.. plot:: _gallery/plot_4_1_4_margenau_hill.py

The Affine class of time-frequency representations
--------------------------------------------------

The Morlet scalogram of two Gaussian atoms.

.. plot:: _gallery/plot_4_2_2_morlet_scalogram_atoms.py

The reassignment method
-----------------------

Reassigned spectrogram of a signal containing sinusoidal and hyperbolic
frequency modulation.

.. plot:: _gallery/plot_4_3_2_reassigned_spectrogram.py

Comparison of different time-frequency distributions and their reassigned
versions. Input signal contains linear, sinusoidal and constant frequency
modulated components.

.. plot:: _gallery/noplot/noplot_4_3_5_wv_smoothed_reassigned.py
.. plot:: _gallery/noplot/noplot_4_3_5_morlet_scalogram_spectrogram.py
.. plot:: _gallery/noplot/noplot_4_3_5_page_margenau_hill.py

Extracting information with the Hough transform
-----------------------------------------------

WVD of a noisy chirp and it's Hough transform. The peak corresponding to the
noisy signal can be see.

.. plot:: _gallery/plot_5_4_2_wv_noisy_chirp.py
.. plot:: _gallery/plot_5_4_2_hough_noisy_chirp.py

The WVD of two simultaneous chirps, followed by it's Hough transform. The WVD
shows a lot of interference terms, but the Hough transform represents this
interference as small side lobes.

.. plot:: _gallery/plot_5_4_2_wv_simultaneous_chirp.py
.. plot:: _gallery/plot_5_4_2_hough_simultaneous_chirp.py
