# CASHEW
[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=symbiotic-engineering/CASHEW)
![GitHub License](https://img.shields.io/github/license/symbiotic-engineering/cashew)

Repository for CASHEW project, submitted for the [DOE Power at Sea Prize](https://www.herox.com/PowerAtSea).

The CASHEW (CArbon Sequestration Harnessing Energy from Waves) system uses wave power to directly pump liquid carbon dioxide into the seabed,
reducing the cost of carbon sequestration. This repository contains the models used as well as many of the figures submitted to the prize.

Read our technical report: [`pubs/CASHEW_power_at_sea_prize.pdf`](https://github.com/symbiotic-engineering/CASHEW/tree/main/pubs/CASHEW_power_at_sea_prize.pdf).

Watch our video: 

[![IMAGE ALT TEXT](https://img.youtube.com/vi/GPyi_EUzb6E/0.jpg)](http://www.youtube.com/watch?v=GPyi_EUzb6E "SEA Lab: Power at Sea CONCEPT Phase Submission - CASHEW")

Run our model:
- To run the Simscape model for transient pressure and flow, simply run the `runCASHEW.m` file in the `src/simscape` folder. Then plot the results using the `plot_results.m` file in the `src/plotting` folder. Note that you must have WEC-Sim installed on your device for this simulation to work. Visit https://wec-sim.github.io/WEC-Sim/dev/user/getting_started.html#download-wec-sim to install. This model is only confirmed to work with MATLAB R2023b.
- To run the MATLAB model for steady state pressure and flow, run the `pressureTempCalc.m` file in the `src/simple_model` folder.
- To run the factor of merit model, run the `fom_sweep.m` file in the `src/econ` folder.

Funding Acknowledgement

This material is based upon work supported by the National Science Foundation Graduate Research Fellowship under Grant No. DGE–2139899. Any opinion, findings, and conclusions or recommendations expressed in this material are those of the authors(s) and do not necessarily reflect the views of the National Science Foundation.
