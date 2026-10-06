## Implement Spin

This work broadly follows Qian (2025) thesis who has implemented spin and polarisation to OSIRIS. 

The code framework borrows from Holger's repo of EPOCH who has partially implemented spin housekeeping and T-BMT into EPOCH, but photon and QED have not been touched.

The code is also based on Chris Arran's work on particle splitting and Yutong's work on linear Breit-Wheeler (the latter is alread on main)

A decision has been made on coding anomalous magnetic dipole moment for different species: for electrons, positrons, muons, and protons, the constant can easily be coded (and is required) in constants.F90. However, for non-proton ion species (e.g. Carbon), an individualised constant is required, and there isn't a general formula to calculate such. Coding a table can be very time consuming and eventually it will require the user to specify the name of the ion species in input.deck. So, instead, I've decided to enable input.deck species block anomalous magnetic dipoel moment parsing for these ion species. THe electron, positrin, muon, and proton ones are hard-coded. If no data is supplied, spin evolution of that species is disabled. Another issue is that to retain unitarity, the electron anomalous magnetic dipole moment should be chi-dependent rather than a constant. At the early stage of code development, I have left that as a constant, but later I will code the chi-dependent a_e(chi_e).

T-BMT: I think Holger has one coefficient wrong (or Qian). The T-BMT I implemented is based on Qian (2025) (2.45) and (2.46) .Importantly, the spin is in the rest frame of the particles, while the fields are in the lab frame. (fac????)

T-BMT is implemented on 6 Oct 2026. Haven't done tests yet. 