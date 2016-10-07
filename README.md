# Residue-Geometry-Network
Amino acid networks (AANs) abstract the protein structure by recording the amino acid contacts and can provide insight into protein function. The Residue Geometry Network (RGN) is a novel AAN construction technique that employs the rigidity analysis tool, FIRST (http://azte.technologypublisher.com/technology/7041), to build the AAN. We have shown that this new construction can be combined with network theory methods  to include the effects of allowed conformal motions and local chemical environments. Importantly, this is done without costly molecular dynamics simulations required by other AAN-related methods, which allows us to analyse large proteins and/or data sets. Random walk simulations using the RGN were also successful in identifying allosteric residues in proteins involved in GPCR signalling (http://www.nature.com/articles/srep33213). 

# Files
RGN.in - text script containing key words to be passed to the main script  
RGN.py - main python script for building RGN and running network analysis

# Required Programs and Packages

Mathematica - https://www.wolfram.com/mathematica/trial/  
NetworkX - https://networkx.github.io/  
FIRST -  http://flexweb.asu.edu   
Reduce - http://molprobity.biochem.duke.edu/index.php?MolProbSID=n45i9vcsk92jg1s30v68fjc5e4&eventID=58  
Reduce Dictionary - http://kinemage.biochem.duke.edu/software/reduce.php  

# KEY WORDS #

provide_chain   
  if n  the chain will be provided in the file 'pdb_chain_list' with format pdb_name chain e.g. 1RRY A for chain A of 1RRY  
  or if y whether to provide the chain in the space below provide_chain using variable chain. This allows the user to define multiple chains  
  
RESTART  
  if RGN has already been run this will take the already formated file and begin again  
provided  
  if y then pdb file will be provided, note there must be a EXPDTA line and REPRESENTATIVE ENSEMBLE line  
  if n then pdb file will be downloaded  
provide_add_chain  
  some pdb files do not have a chain column, if y this will ad one termed A  
h_added   
  if y then reduce won't add hydrogens  
  note most pdb files do not have hydrogens  
water  
  name of water in file  
het_remove  
  list of heteroatoms to be removed from analysis, these are often crystal artifacts  
dilution_plot   
  if y then calculate dilution plot (note this will take longer than normal FIRST analysis)  
RC_identify   
  if y then calculate rigid clusters are a hydron bond energy cutoff of %RC_cutoff  
RC_cutoff   
  see above  
RC_atom  
  atoms used for rigid clusters  

phobe_color   
  hydrophobic interactions color in graph  
cov_color  
   covalent interactions color in graph   
hb_color   
  hydrogen bond interaction color in graph  
phobe_strength   
  hydrophobic strength in weighted network, note this was identified in "Residue Geometry Networks" paper in Scientific Reports  
cov_strength  
  hydrophobic strength in weighted network, note this was identified in "Residue Geometry Networks" paper in Scientific Reports  

deg_cutoff  
  Hcut when building graph to  calculate degree centrality  
bet_cutoff    
  Hcut when building graph to  calculate betweenness centrality    
clo_cutoff  
  Hcut when building graph to  calculate closeness centrality    
no_phobe   
  if y then phobe interactions will be removed from  closeness centrality measurement and matrix if matrix == 'y'   

matrix_set   
  if y then run EVT analysis  
EVT_centrality   
  if y then calculate centarlity of EVT network    
dist_matrix_restart   
  if distance matrix has already been calculated then set y and avoid recalculating it    
evt_matrix_restart   
  if EVT matrix has already been calculated then set y and avoid recalculating it  
Nmax= 150  
mat_cutoff    
  matrix cutoff for EVT analysis  

codons  
  run codon analysis  
codons_path 
  path to codons list  
  For this analysis biopython paml must be downloaded  
  and a codons list must be supplied  
  FORMAT  
       1i41.A.YAL012W.core     1       Y       ATG:ATG:ATG:ATG  
       1i41.A.YAL012W.core     2       A       ATG:ATG:ATG:ATG  


shortest_path  
  calculate shortest between residues in file f=open('%s_shortest_path.txt'%pdb)  
  format:  
  begin A 1 THR  
  end A 4 THR  
makeplot   
  plot degree against closeness centrality in file clo_deg_scatter.png  
subgraph   
  analyse subgraph  
sublist   
  which nodes in subgraph_list  

net_vis   
  visualise network  
node_col_opt   
  if y then color nodes in node_col list a dif colour  
node_col_list   
  residues to be coloured a dif colour   
node_alt_col   
  alternative colouring for residues  


# ADDING INTERACTIONS MANUALLY 
It is imporatant that the interactions be inspected, in particular when crystal  structures of low resolution are used.  
This can lead to an abnormally short distance between the interacting groups, for which the energy function used in FIRST will penalise the interaction, assigning a positive energy to the salt bridge. In these cases, the hbond.out can be modified, and modest interaction energies can be assigned. Then, move the hbond.out file to hbond.in ( as well as hphobe.out to hphobe.in and cov.out to cov.in) and switch the keyword argument write_interactions to read_interactions in pgn.in 
