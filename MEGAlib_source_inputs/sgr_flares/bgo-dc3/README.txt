2024/June/23 Che-Yen Jimmy Chu
This is the description to the source model of magnetar bursts (SRG J1935+2154 and 1E 1048.1-5937)
2 light curves and 2 different flux level are included

Burst 1: SGR J1935+2154, a bright burst with complex light curve (source file: MgtBurst_bright_complex.source)
Burst 2: 1E 1048.1-5937, a dim burst with simple light curve (source file: MgtBurst_dim_simple.source)
Burst 3: spectrum of SGR J1935 burst + light curve of 1E 1048 burst (source file: MgtBurst_bright_simple.source)
==> a bright burst with simple light curve
Burst 4: spectrum of 1E 1048 burst + light curve of SGR J1935 burst (source file: MgtBurst_dim_complex.source)
==> a dim burst with complex light curve

Reference:
SGR J1935+2154: Li, et al., 2021, NatAstron., 5, 378
A bright burst with radio conterpart (FRB) confirmed
Light curve: using 10-30 keV light curve since 27-250 keV light curve suffers pile-up
Spectrum: using cutoff power-law with power index 1.56 and cutoff energy 83.89

1E 1048.1-5937: Gavriil, et al., 2002, Nature, 419, 142
The first burst which was discovery from an AXP (not from SGR)
Light curve: using 2-20 keV light curve
Spectrum: using cutoff power-law with power index 0.89 and assumed* cutoff energy 37.01
*The spectrum fitting reported is up to 40 keV, which makes the cutoff hard to idnetify, so I choose the avereged cutoff energy of magnetar burst (Collazzi, et al., 2015, ApJS, 218, 11).

Note:
These models with a lower energy range (10-500 keV) are for future data challenges which we may use shields event for analysis. For DC3, please refer to the previous directory for the models with the standard energy range (0.1-10 MeV).
