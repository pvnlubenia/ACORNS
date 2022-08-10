=================================================================================

  ACORNS: Absolute COncentration Robustness of Networks of Shinar-feinberg type

=================================================================================

Matlab was used to develop the functions used here.



=========
Functions
=========The function acr returns a list of species (together with the Shinar-Feinberg (SF) pair associated with each and the deficiency of the building block subnetwork containing the SF-pair) with absolute concentration robustness (ACR) in a chemical reaction network (CRN), if they exist. ACR in a species is checked for each SF-pair even if the species is already determined to have ACR considering a different SF-pair. If no species is found or the network is not of SF-type, a message appears saying so. Furthermore, the output variables 'model', 'R', 'F', and 'ACR_species' allow the user to view the following, respectively:

   - Complete network with all the species listed in the 'species' field of the structure 'model'
   - Matrix of reaction vectors of the network
   - Kinetic order matrix of the network
   - List of species with ACR

acr uses the following functions:     1. model_species
           - OUTPUT: Creates a complete model structure where the species list is filled out. The output variables 'model' and 'm' allow the user to view the complete model and the number of species, respectively.
           - INPUT: model: a structure representing the CRN (see details below)

     2. stoich_matrix           - OUTPUT: Creates the stoichiometric matrix of the system. The output variables 'N', 'reactant_complex', 'product_complex', and 'r' allow the user to view the stoichiometric matrix, matrix of reactant complexes, matrix of product complexes, and number of reactions, respectively.
           - INPUTS
                - model: a structure representing the CRN (see details below)
                - m: number of species
     3. kin_ord_matrix           - OUTPUT: Returns the kinetic order matrix of a system. It is assumed that the system has mass action kinetics.
           - INPUTS
                - model: a structure representing the CRN (see details below)
                - m: number of species
     4. R_is_span_union
           - OUTPUT: Returns a value of 1 for the variable binary_decomp if R is in the union of span(B1) and span(B2), i.e., an independent binary decomposition of R is formed. The output variable 'span_B1' allows the user to view the elements (reaction numbers) in span(B1).
           - INPUTS
                - B1: an array of the reaction numbers of an SF-pair
                - B2: an array of the reaction numbers of the vectors added to extend B1 to a basis for R                - R: a matrix of reaction vectors of a CRN
     5. deficiency
           - OUTPUT: Creates a network out of a list of reactions from another network, then returns the deficiency value of the newly-formed CRN. The output variables 'model_N1' and delta1 allow the user to view the new network formed and the deficiency of the new network, respectively.
           - INPUTS
                - model: a structure representing the CRN (see details below)
                - span_B1: list of reaction numbers of 'model' that will be used to form model_N1

     6. is_PL_RDK
           - OUTPUT: Returns a value of 1 if the power law system has reactant-determined kinetics (PL-RDK), 0 otherwise.
           - INPUTS
                - model: a structure representing the CRN (see details below; it is assumed that the CRN has power law kinetics and that the field 'species' have already been filled out by another function)
                - m: number of species

     7. is_weakly_reversible
           - OUTPUT: Returns a value of 1 if the network is weakly reversible, 0 otherwise.
           - INPUTS
                - model: a structure representing the CRN (see details below; it is assumed that the CRN has power law kinetics and that the field 'species' have already been filled out by another function)
                - m: number of species

     8. add_vertices
           - OUTPUT: Adds all the complexes as vertices in a (di)graph. The output variables 'g' and 'vertices' allow the user to view the (di)graph with vertices and the list of vertices, respectively.
           - INPUTS
                - g: (di)graph without vertices
                - n: number of complexes
                - Y: matrix of complexes
                - model: a structure representing the CRN (see details below; it is assumed that the CRN has power law kinetics and that the field 'species' have already been filled out by another function)

     9. add_edges           - OUTPUT: Adds all the edges in a (di)graph. The output variables 'g' and 'edges' allow the user to view the (di)graph with edges and the list of edges, respectively.
           - INPUTS
                - g: (di)graph without edges
                - vertices: list of vertices
                - r: number of reactions
                - Y: matrix of complexes
                - reactant_complex: matrix of reactant complexes
                - product_complex: matrix of product complexes

     10. in_same_linkage_class
           - OUTPUT: Returns a value of 1 if the complexes are in the same linkage class in 'model', 0 otherwise.
           - INPUTS
                - complex1, complex2: strings representing the complexes we want to check
                - model: a structure representing the CRN (see details below; it is assumed that the field 'species' have already been filled out by another function)
                - m: number of species
     11. is_nonterminal
           - OUTPUT: Returns a value of 1 if the complex is nonterminal in 'model', 0 otherwise.
           - INPUTS
                - check_complex: string representing the complex we want to check
                - model: a structure representing the CRN (see details below; it is assumed that the field 'species' have already been filled out by another function)
                - m: number of species


====
Note
====

acr2 is the same as acr but returns each species (together with the SF-pair associated with it and the deficiency of the building block block subnetwork containing the SF-pair) with ACR as they are found, if they exist. The time elapsed is also shown per species found and the end of the running time. This is perfect for networks with high rank and a lot of SF-pairs. An alphabetical list of the species is also returned at the end.



=================================
How to fill out 'model' structure
=================================

'model' is the input for the function acr. It is a structure, representing the CRN, with the following fields:

   - id: name of the model
   - species: a list of all species in the network; this is left blank since incorporated into the function is a step which compiles all species used in the model
   - reaction: a list of all reactions in the network, each with the following subfields:
        - id: a string representing the reaction
        - reactant: has the following further subfields:
             - species: a list of strings representing the species in the reactant complex
             - stoichiometry: a list of numbers representing the stoichiometric coefficient of each species in the reactant complex (listed in the same order of the species)
        - product: has the following further subfields:
             - species: a list of strings representing the species in the product complex
             - stoichiometry: a list of numbers representing the stoichiometric coefficient of each species in the product complex (listed in the same order of the species)
        - reversible: has the value true or false indicating if the reaction is reversible or not, respectively
        - kinetic: has the following further subfields:
             - reactant1: a list of numbers representing the kinetic order of each species in the reactant complex in the left to right direction (listed in the same order of the species)
             - reactant2: a list of numbers representing the kinetic order of each species in the reactant complex in the right to left direction (listed in the same order of the species) (empty if the reaction is not reversible)

To fill out the 'model' structure, write a string for 'model.id': this is just to put a name to the network. To add the reactions to the network, use the function addReaction where the output is 'model'. addReaction is developed to make the input of reactions of the CRN easier than the input in [7]:

   addReaction
      - OUTPUT: Returns a structure called 'model' with added field 'reaction' with subfields 'id', 'reactant', 'product', 'reversible', and 'kinetic'. The output variable 'model' allows the user to view the network with the added reaction.
      - INPUTS
           - model: a structure, representing the CRN
           - id: visual representation of the reaction, e.g., reactant -> product (string)
           - reactant_species: species of the reactant complex (cell)
           - reactant_stoichiometry: stoichiometry of the species of the reactant complex (cell)
           - reactant_kinetic: kinetic orders of the species of the reactant complex (array)
           - product_species: species of the product complex (cell)
           - product_stoichiometry: stoichiometry of the species of the product complex (cell)
           - product_kinetic: "kinetic orders" of the species of the product complex, if the reaction is reversible (array); if the reaction in NOT reversible, leave blank
           - reversible: logical; whether the reaction is reversible or not (true or false)
      * Make sure the function addReaction is in the same folder/path being used as the current working directory.

Note that for the functions acr and acr2:

     1. It is assumed that the CRN has a positive equilibrium.
     2. It is also assumed that the CRN has power law kinetics.
     3. The CRN should have at least 2 species and 2 reactions (to form an SF-pair).
     4. Notes 2 and 3 imply that we assume the CRN is a power law system of SF-type.
     5. This code is based largely on [2] with modifications based on [5].



========
Examples
========

8 examples are included in this folder:

   - Example 1: Early STAT signaling network coupled with the receptor complex formation upon interferon induction [2]

   - Example 2: A network with 6 components [1]

   - Example 3: A simple two-species mass-action reaction system [4]

   - Example 4: Power law approximation of the pre-industrial carbon cycle model [3]

   - Example 5: A network with only one species [6]

   - Example 6: Another network with only one species [6]

   - Example 7: A network with two species [6]

   - Example 8: A network with power law kinetics [3]



===================
Contact Information
===================

For questions, comments, and suggestions, feel free to contact me at pvnlubenia@yahoo.co.uk.


- Patrick Lubenia (18 July 2022)



==========
References
==========

   [1] Eloundou-Mbebi J, KŸken A, Omranian N, Kleesen S, Neigenfind J, Basler G, Nikoloski Z (2016) A network property necessary for concentration robustness. Nat Commun 7(13255):1-7. https://doi.org/10.1038/ncomms13255

   [2] Fontanil L, Mendoza E, Fortun N (2021) A computational approach to concentration robustness in power law kinetic systems of Shinar-Feinberg type. MATCH Commun Math Comput Chem 86(3):489-516.

   [3] Fortun N, Mendoza E (2020) Absolute concentration robustness in power law kinetic systems. MATCH Commun Math Comput Chem 85(3):669-691.

   [4] Kuwahara H, Umarov R, Almasri I, Gao X (2017) ACRE: Absolute concentration robustness exploration in module-based combinatorial networks. Synth Biol 2(1):1-6. https://doi.org/10.1093/synbio/ysk01

   [5] Lao A, Lubenia P, Magpantay D, Mendoza E (2022) Concentration robustness in LP kinetic systems. MATCH Commun Math Comput Chem 88(1):29-66. https://doi.org/10.46793/match.88-1.029L

   [6] Meshkat N, Shiu A, Torres A (2021) Absolute concentration robustness in networks with low-dimensional stoichiometric subspace. Vietnam J Math https://doi.org/10.1007/s10013-021-00524-5

   [7] Soranzo N, Altafini C (2009) ERNEST: a toolbox for chemical reaction network theory. Bioinform 25(21):2853-2854. https://doi.org/10.1093/bioinformatics/btp513