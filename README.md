## fbdR is an R package for estimating speciation and extinction rates from phylogenetic trees and fossil occurrence data

### Estimating rates from extant phylogenies and stratigraphic ranges

Available methods include

* birth-death process model (Stadler, 2012, eqn. 2) 
* birth-death process model (Kieding, 1975; Silvestro et al. 2014, eqn. 9)

Described here

* Stadler, T. 2012. How can we improve accuracy of macroevolutionary rate estimates? Syst Biol, 62: 321–329, 2012.
* Keiding N. 1975. Maximum likelihood estimation in the birth-death process. Ann. Stat. 3:363–372.
* Silvestro, D et al. 2014. Bayesian Estimation of Speciation and Extinction from Incomplete Fossil Occurrence Data. Syst Bio 63: 349-367.

The implementation of these methods available in the `fbdR` package assume that the age of all speciation and extinction events are known. Estimates of speciation and extinction can be obtained using maximum likelihood for trees and stratigraphic ranges, independently or jointly.

### Estimating rates from fossil occurrence data

Available methods include

* boundary crosser method (Foote, 2000)
* three-timer method (Alroy 2008)
* gap-filler method (Alroy, 2014)

Desribed here 

* Foote, M. 2000. Origination and extinction components of taxonomic diversity: General problems. Paleobiology 26: 74-102.
* Alroy, J et al. 2008. Dynamics of origination and extinction in the marine fossil record. PNAS 105: 11536-11542. 
* Alroy, J 2014. Accurate and precise estimates of origination and extinction rates. Paleobiology 40: 374-397.

### Package installaton

The latest version can be installed in R using the package devtools:

    library(devtools)
    install_github("rachelwarnock/fbdr")

`fbdR` was designed to work with output from the simulation packages `TreeSim` and `FossilSim` but can be used to estiamte diversification rates using any simulated or empirical data, so long as they match the input format.
    
### For further information and examples see the package documentation