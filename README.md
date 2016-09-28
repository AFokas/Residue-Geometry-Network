# Residue-Geometry-Network
Amino acid networks (AANs) abstract the protein structure by recording the amino acid contacts and can provide insight into protein function. The Residue Geometry Network (RGN) is a novel AAN construction technique that employs the rigidity analysis tool, FIRST (http://azte.technologypublisher.com/technology/7041), to build the AAN. We have shown that this new construction can be combined with network theory methods  to include the effects of allowed conformal motions and local chemical environments. Importantly, this is done without costly molecular dynamics simulations required by other AAN-related methods, which allows us to analyse large proteins and/or data sets. Random walk simulations using the RGN were also successful in identifying allosteric residues in proteins involved in GPCR signalling (http://www.nature.com/articles/srep33213). 

RGN.in - text script containing key words to be passed to the main script.
RGN.py - main python script for building RGN and running network analysis

Required Programs and Packages:

Mathematica - https://www.wolfram.com/mathematica/trial/
NetworkX - https://networkx.github.io/
FIRST -  http://flexweb.asu.edu 
Reduce - http://molprobity.biochem.duke.edu/index.php?MolProbSID=n45i9vcsk92jg1s30v68fjc5e4&eventID=58
Reduce Dictionary - http://kinemage.biochem.duke.edu/software/reduce.php

##### ADDING INTERACTIONS MANUALLY ####
It is imporatant that the interactions be inspected, in particular when crystal  structures of low resolution are used.  
This can lead to an abnormally short distance between the interacting groups, for which the energy function used in FIRST will penalise the interaction, assigning a positive energy to the salt bridge. In these cases, the hbond.out can be modified, and modest interaction energies can be assigned. Then, move the hbond.out file to hbond.in ( as well as hphobe.out to hphobe.in and cov.out to cov.in) and switch the keyword argument write_interactions to read_interactions in pgn.in 
