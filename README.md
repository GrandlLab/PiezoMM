# PiezoMM

[![Build Status](https://ci.appveyor.com/api/projects/status/github/neuro-myoung/PiezoMM.jl?svg=true)](https://ci.appveyor.com/project/neuro-myoung/PiezoMM-jl)
[![Coverage](https://codecov.io/gh/neuro-myoung/PiezoMM.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/neuro-myoung/PiezoMM.jl)
[![Coverage](https://coveralls.io/repos/github/neuro-myoung/PiezoMM.jl/badge.svg?branch=main)](https://coveralls.io/github/neuro-myoung/PiezoMM.jl?branch=main)

## Description
This is a julia package written to simulate Piezo1 channel activity in response to an *in silico* tension clamp. 

Briefly, the channels are randomly seeded onto a grid representing the cell membrane with a probability such that channel densities approximate those determined in Lewis 2021 (~100 channels/um^2). The image below shows an example of channel seeding for a single simulation.

![channels](/assets/channels.png)

The tension stimulus is applied and allowed to diffuse following a normalized Gaussian distribution. Tension diffusion rates are bound by the findings of Shi 2018 and Shi 2021. 

![tension](/assets/tension.gif) kinetic_model

Each individual channel's behavior is modelled as a discrete Markov process. The transition probabilities are derived from the 4-state kinetic model of Piezo1 determined in Lewis 2017.

![rate_model](/assets/kinetic_model.png)

An example simulation can be seen below:
![rate_model](/assets/ex_sim.png)

A more detailed interactive example can be found in the Pluto Notebook included in the repository.

## References
1. Lewis AH, Cui AF, McDonald MF, Grandl J. 2017. Transduction of repetitive mechanical stimuli by piezo1 and piezo2 ion channels. Cell Rep 19:2572–2585. doi:10.1016/j.celrep.2017.05.079

2. Lewis AH, Grandl J. 2021. Piezo1 ion channels inherently function as independent mechanotransducers. eLife 10. doi:10.7554/eLife.70988

3. Shi Z, Graber ZT, Baumgart T, Stone HA, Cohen AE. 2018. Cell membranes resist flow. Cell 175:1769-1779.e13. doi:10.1016/j.cell.2018.09.054

4. Shi Z, Innes-Gold S, Cohen AE. 2022. Membrane tension propagation couples axon growth and collateral branching. BioRxiv. doi:10.1101/2022.01.09.475560

## Prerequisites

Before you begin make sure you have Julia 1.7 or higher. It may work with other versions but I have not tested it.

## Installation

```
git clone https://github.com/neuro-myoung/PiezoMM.git
```

All the required packages are contained in the Manifest. Next, the following two commands will activate the environment and install the requisite packages. Make sure you type these commands from the package manager interface that can be accessed through julia by pressing `]` in the command line.

```
(v1.0) pkg> activate
(PiezoMM) pkg> instantiate
``` 

## Contributing
To contribute to **PiezoMM**, follow these steps:

1. Fork this repository.
2. Create a branch: git checkout -b *branch_name*.
3. Make your changes and commit them: git commit -m '*commit_message*'
4. Push to the original branch: git push origin *project_name* *location*
5. Create a pull request.

Alternatively see the GitHub [documentation](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request) on creating a pull request.

## Contributors

[@neuro-myoung](https://github.com/neuro-myoung)

[@msindoni] (https://github.com/msindoni)

## Contact

If you want to contact me you can reach me at michael.young@duke.edu

## License
This project uses an [MIT License](https://opensource.org/licenses/MIT)

