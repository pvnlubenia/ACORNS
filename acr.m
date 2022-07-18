% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                           %
%    acr                                                                    %
%                                                                           %
%                                                                           %
% OUTPUT: Returns a list of species (together with the Shinar-Feinberg (SF) %
%    pair associated with each and the deficiency of the building block     %
%    subnetwork containing the SF-pair) with absolute concentration         %
%    robustness (ACR) in a chemical reaction network (CRN), if they exist.  %
%    ACR in a species is checked for each SF-pair even if the species is    %
%    already determined to have ACR considering a different SF-pair. If no  %
%    species is found or the network is not of SF-type, a message appears   %
%    saying so. The output variables 'model', 'R', 'F', and 'ACR_species'   %
%    allow the user to view the following, respectively:                    %
%       - Complete network with all the species listed in the 'species'     %
%            field of the structure 'model'                                 %
%       - Matrix of reaction vectors of the network                         %
%       - Kinetic order matrix of the network                               %
%       - List of species with absolute concentration robustness            %
%                                                                           %
% INPUT: model: a structure, representing the CRN (see README.txt for       %
%    details on how to fill out the structure)                              %
%                                                                           %
% Notes:                                                                    %
%    1. It is assumed that the CRN has a positive equilibrium.              %
%    2. It is also assumed that the CRN has power law kinetics.             %
%    3. The CRN should have at least 2 species and 2 reactions (to form an  %
%          SF-pair).                                                        %
%    4. Notes 2 and 3 imply that we assume the CRN is a power law kinetic   %
%          system of SF-type.                                               %
%    5. This code is based largely on [1] with modifications based on [2].  %
%                                                                           %
% References                                                                %
%    [1] Fontanil L, Mendoza E, Fortun N (2021) A computational approach to %
%           concentration robustness in power law kinetic systems of        %
%           Shinar-Feinberg type. MATCH Commun Math Comput Chem             %
%           86(3):489-516.                                                  %
%    [2] Lao A, Lubenia P, Magpantay D, Mendoza E (2022) Concentration      %
%           robustness in LP kinetic systems. MATCH Commun Math Comput Chem %
%           88(1):29-66. https://doi.org/10.46793/match.88-1.029L           %
%    [3] Soranzo N, Altafini C (2009) ERNEST: a toolbox for chemical        %
%           reaction network theory. Bioinform 25(21):2853–2854.            %
%           https://doi.org/10.1093/bioinformatics/btp513                   %
%                                                                           %
% Created: 17 June 2022                                                     %
% Last Modified: 18 July 2022                                               %
%                                                                           %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



function [model, R, F, ACR_species] = acr(model)
    
    %
    % Step 1: Create a list of all species indicated in the reactions
    %
    
    [model, m] = model_species(model);
    
    
    
    %
    % Step 2: Form stoichiometric matrix N
    %
    
    [N, reactant_complex, ~, r] = stoich_matrix(model, m);
    
    
    
    %
    % Step 3: Form kinetic order matrix F
    %    Reminder: Algorithm assumes the system is mass action
    %
    
    F = kin_ord_matrix(model, m);
    
    
    
    %
    % Step 4: Get the matrix of reaction vectors of the network and its rank
    %
    
    R = N';
    
    
    
    %
    % Step 5: Create a list of reactant complexes
    %
    
    % Initialize list of reactant complexes
    reactant_complex_list = { };
    
    % Go through each reactant complex
    for i = 1:size(reactant_complex, 2)
        
        % Get each column of reactant_complex
        reactant_complex_ = reactant_complex(:, i);

        % Get the index numbers of nonzero entries
        reactant_complex_nnz = find(reactant_complex_);
        
        % For the zero complex
        if numel(reactant_complex_nnz) == 0
            complex = '0';
        
        % Otherwise
        else
            
            % Check which species appear in the complex
            for j = 1:numel(reactant_complex_nnz)
                
                % For the first species
                if j == 1
                    
                    % Don't show the stoichiometry if it's 1
                    if reactant_complex_(reactant_complex_nnz(j)) == 1
                        complex = [model.species{reactant_complex_nnz(j)}];
                    
                    % Otherwise, include it
                    else
                        complex = [num2str(reactant_complex_(reactant_complex_nnz(j))), model.species{reactant_complex_nnz(j)}];
                    end
                
                % We need the + sign for succeeding species
                else
                    if reactant_complex_(reactant_complex_nnz(j)) == 1
                        complex = [complex, '+', model.species{reactant_complex_nnz(j)}];
                    else
                        complex = [complex, '+', num2str(reactant_complex_(reactant_complex_nnz(j))), model.species{reactant_complex_nnz(j)}];
                    end
                end
            end 
        end
        
        % Add the complex in the list
        reactant_complex_list{end+1} = complex;
    end
    
    
    
    %
    % Step 6: Get all Shinar-Feinberg pairs
    %
    
    % Initialize list of Shinar-Feinberg pairs
    SF_pair = [ ];
    SF_pair_id = { };
    
    % Go through each kinetic order vector
    for i = 1:r
        
        % Compare to each of the other kinetic order vectors
        for j = i+1:r
            
            % Get the absolute values of all entries for each row, then add the rows
            nonzeros = abs(F(i, :)) + abs(F(j, :));
            
            % Make all nonzero entries equal to 1 (for comparison later with true or false values)
            nonzeros(nonzeros ~= 0) = 1;
            
            % Check which kinetic orders match
            kinetics = F(i, :) == F(j, :);
            
            % Ignore if both kinetic orders are 0
            for k = 1:m
                if (F(i, k) == 0 && F(j, k) == 0)
                    kinetics(k) = 0;
                end
            end
            
            % Get only those rows that differ by 1 kinetic order
            if nnz(nonzeros) - 1 == sum(kinetics)
                
                % Get the index of the kinetic order that is different
                index = find(~(nonzeros == kinetics));
                
                % Compile Shinar-Feinberg pairs in matrix form: reactions i and j form a Shinar-Feinberg pair in species 'index')
                SF_pair(end+1, :) = [i, j, index];
                
                % Get the name of the species corresponding to that kinetic order
                species = model.species{index};
                
                % Compile a list of Shinar-Feinberg pairs
                SF_pair_id{end+1} = ['(R' num2str(i) ', R' num2str(j) ') in ' species];
            end
        end
    end
    
    % If there are no Shinar-Feinberg pairs, exit the algorithm
    if size(SF_pair, 1) == 0
        ACR_species = { };
        fprintf([model.id ' is not of Shinar-Feinberg type. The algorithm cannot be used. \n\n']);
        return
    end
    
    
    
    %
    % Step 7: Form a basis for the rowspace of R
    %
    
    % Write R in reduced row echelon form: the transpose of R is used so basis_reac_num will give the pivot rows of R
    %    - basis_reac_num gives the row numbers of R which form a basis for the rowspace of R
    [~, basis_reac_num] = rref(R');
    
    % Form the basis
    basis = R(basis_reac_num, :);
    
    % Initialize list of species with absolute concentration robustness (ACR)
    % This also serves as a control if no independent binary decomposition is found: it will remain empty
    ACR_species = { };
    
    
    
    %
    % Step 8: Get a Shinar-Feinberg pair
    %
    
    % Initialize flag for the loops
    flag = 0;
    
    % Go through each Shinar-Feinberg pair
    for i = 1:size(SF_pair, 1)
        
        % Get the reaction numbers of the Shinar-Feinberg pair
        SF_pair1 = SF_pair(i, 1);
        SF_pair2 = SF_pair(i, 2);
        
        
        
    %
    % Step 9: Extend the pair to a basis for the rowspace of R
    %
        
        % Note: A vector v is a linear combination of vectors in a matrix A if the rank of the matrix [A, v] (v is appended to A) is the same as the rank of A
        
        % Initialize the extended basis
        basis_SF = [ ];
        
        % Case 1: The Shinar-Feinberg pair is NOT linearly independent, i.e., each is a linear combo of the other
        if rank([R(SF_pair1, :)', R(SF_pair2, :)']) == rank(R(SF_pair1, :)')
        
            % Use the first reaction vector of the Shinar-Feinberg pair to start the formation of the basis
            basis_SF(end+1, :) = R(SF_pair1, :);
            
            % This forms B1
            B1 = SF_pair1;
            
        % Case 2: The Shinar-Feinberg pair is linearly independent
        else
        
            % Use the reaction vectors to start the formation of the basis
            basis_SF(end+1, :) = R(SF_pair1, :);
            basis_SF(end+1, :) = R(SF_pair2, :);
            
            % These form B1
            B1 = [SF_pair1, SF_pair2];
        end
        
        % Initialize list of reaction vectors added to extend basis_SF to a basis for R
        added_reaction = [ ];
        
        % Keep on going through each vector in 'basis' until there are s [= rank(N)] elements in basis_SF
        while size(basis_SF, 1) < size(basis, 1)
        
            % Go through each vector in 'basis'
            for j = 1:size(basis, 1)
            
                % If the reaction vector is NOT a linear combination of basis_SF
                if rank([basis_SF', basis(j, :)']) ~= rank(basis_SF')
                
                    % Add this vector to basis_SF
                    basis_SF(end+1, :) = R(basis_reac_num(j), :);
                    
                    % Take note which reaction vectors from 'basis' are added to basis_SF
                    added_reaction(end+1) = j;
                    
                    % Stop forming basis_SF when its size reaches the rank s of the network
                    % Note: basis_SF will also have a rank s since it is composed of linearly independent vectors that span R
                    if size(basis_SF, 1) == size(basis, 1)
                        break
                    end
                end
            end
        end    
        
        % Form the second set of the separated basis vectors
        B2 = basis_reac_num(added_reaction);
        
        
        
    %
    % Step 10: Check if R is the union of span(B1) and span(B2)
    %
        
        % Form span_B1 and span_B2
        % binary_decomp is 1 if R is the union of span(B1) and span(B2), 0 otherwise
        [binary_decomp, span_B1] = R_is_span_union(B1, B2, R);
        
        % Initialize list of power sets of B2 (for use in case we don't get an independent binary decomposition)
        power_set_B2 = { };
        
        % Create a list of power sets of B2
        % We exclude the empty set and the full set
        for j = 1:numel(B2)-1
            power_set_B2{j} = nchoosek(B2, j);
        end
        
        % Case 1: R is NOT the union of span(B1) and span(B2)
        if binary_decomp == 0
            
            % NOTE: In the interest of understanding the process, I intentionally did not make the succeeding iterated loops into a function even if it will be repeated later
            
            
            
    %
    % Step 11: Transfer some elements of B2 to B1 until an independent binary decomposition is found
    %
            
            % Go through each set in the power set
            for j = 1:numel(power_set_B2)
                
                % Go through each element of a set in the power set
                for k = 1:size(power_set_B2{j}, 1)
                    
                    % Transfer elements from B2 to B1
                    B1_new = [B1 power_set_B2{j}(k, :)];
                    B2_new = B2;
                    B2_new(ismember(B2_new, power_set_B2{j}(k, :))) = [ ];
                    
                    % Check if R is the union of span(B1) and span(B2)
                    [binary_decomp, span_B1] = R_is_span_union(B1_new, B2_new, R);
                    
                    % Once successful in getting an independent binary decomposition
                    if binary_decomp == 1
                        
                        
                        
    %
    % Step 12: Get the deficiency of the network induced by span(B1)
    %
                        
                        % Create network N1 called model_N1 from span(B1) and compute its deficiency delta1
                        [model_N1, delta1] = deficiency(model, span_B1);
                        
                        
                        
    %
    % Step 13: Repeat from Step 11 until the deficiency of the network induced by span(B1) is less than or equal to 1
    %
                        
                        % Case 1: The deficiency is 0
                        if delta1 == 0
                            
                            
                            
    %
    % Step 14: For deficiency 0, check if the induced network is a weakly reversible power law system with reactant-determined kinetics (PL-RDK) with the Shinar-Feinberg pair in the same linkage class in the induced network
    %
                            
                            PL_RDK = is_PL_RDK(model_N1, m);
                            if PL_RDK
                                
                                % Check that model_N1 is weakly reversible
                                weakly_reversible = is_weakly_reversible(model_N1, m);
                                if weakly_reversible
                                    
                                    % Get the reactant complex of each reaction in the Shinar-Feinberg pair
                                    reactant_complex1 = reactant_complex_list{SF_pair1};
                                    reactant_complex2 = reactant_complex_list{SF_pair2};
                                    
                                    % Check if the reactant complexes are in the same linkage class in model_N1
                                    same_linkage_class = in_same_linkage_class(reactant_complex1, reactant_complex2, model_N1, m);
                                    
                                    % If the reactant complexes are in the same linkage class
                                    if same_linkage_class
                                        
                                        % Add the species with its associated Shinar-Feinberg pair to the list ACR_species
                                        ACR_species{end+1} = [model.species{SF_pair(i, 3)} ' [SF-pair (R' num2str(SF_pair(i, 1)) ', R' num2str(SF_pair(i, 2)) ') | deficiency ' num2str(delta1) ']'];
                                        
                                        % This is created so we can exit the iterated loops
                                        flag = 1;
                                        
                                        % This breaks out of the loop for k
                                        break
                                    end
                                end
                            end
                        
                        % Case 2: The deficiency is 1
                        elseif delta1 == 1
                             
                          
                          
    %
    % Step 15: For deficiency 1, check if the induced network is a PL-RDK system AND that the reactant complexes of the Shinar-Feinberg pairs are nonterminal complexes
    %
    
                            PL_RDK = is_PL_RDK(model_N1, m);
                            if PL_RDK
                                
                                % Get the reactant complex of each reaction in the Shinar-Feinberg pair
                                reactant_complex1 = reactant_complex_list{SF_pair1};
                                reactant_complex2 = reactant_complex_list{SF_pair2};
                                
                                % Check that the reactant complexes are nonterminal
                                if (is_nonterminal(reactant_complex1, model_N1, m) && is_nonterminal(reactant_complex2, model_N1, m))
                                
                                    % Add the species with its associated Shinar-Feinberg pair to the list ACR_species
                                    ACR_species{end+1} = [model.species{SF_pair(i, 3)} ' [SF-pair (R' num2str(SF_pair(i, 1)) ', R' num2str(SF_pair(i, 2)) ') | deficiency ' num2str(delta1) ']'];
                                    
                                    % This is created so we can exit the iterated loops
                                    flag = 1;
                                    
                                    % This breaks out of the loop for k
                                    break
                                end
                            end
                        end
                    end
                end
                
                % This gets the signal from the loop for k
                if flag == 1
                   
                   % Reset flag
                   flag = 0;
                    
                   % This breaks out of the loop for j so we can move to the next Shinar-Feinberg pair
                   break
                end
            end
      
        % Case 2: R is the union of span(B1) and span(B2)
        % We jump to Step 12 and continue up to Step 15
        else
            
            % Step 12
            [model_N1, delta1] = deficiency(model, span_B1);
            
            % Step 13
            % Case 1: The deficiency is 0
            if delta1 == 0
                
                % Step 14
                PL_RDK = is_PL_RDK(model_N1, m);
                if PL_RDK
                    weakly_reversible = is_weakly_reversible(model_N1, m);
                    if weakly_reversible
                        reactant_complex1 = reactant_complex_list{SF_pair1};
                        reactant_complex2 = reactant_complex_list{SF_pair2};
                        same_linkage_class = in_same_linkage_class(reactant_complex1, reactant_complex2, model_N1, m);
                        if same_linkage_class
                            ACR_species{end+1} = [model.species{SF_pair(i, 3)} ' [SF-pair (R' num2str(SF_pair(i, 1)) ', R' num2str(SF_pair(i, 2)) ') | deficiency ' num2str(delta1) ']'];
                            % No break here since we don't want to exit checking the other pairs
                        end
                    end
                end

            % Case 2: The deficiency is 1
            elseif delta1 == 1
            
                % Step 15
                PL_RDK = is_PL_RDK(model_N1, m);
                if PL_RDK
                    reactant_complex1 = reactant_complex_list{SF_pair1};
                    reactant_complex2 = reactant_complex_list{SF_pair2};
                    if (is_nonterminal(reactant_complex1, model_N1, m) && is_nonterminal(reactant_complex2, model_N1, m))
                        ACR_species{end+1} = [model.species{SF_pair(i, 3)} ' [SF-pair (R' num2str(SF_pair(i, 1)) ', R' num2str(SF_pair(i, 2)) ') | deficiency ' num2str(delta1) ']'];
                        % No break here since we don't want to exit checking the other pairs
                    end
                end
            
            % Case 3: The deficiency is greater than 1
            else
                
                % NOTE: This repeats STEPS 11-15 where we go through each set in the power set until we get a subnetwork with deficiency less than or equal to 1
                % Step 11
                for j = 1:numel(power_set_B2)
                    for k = 1:size(power_set_B2{j}, 1)
                        B1_new = [B1 power_set_B2{j}(k, :)];
                        B2_new = B2;
                        B2_new(ismember(B2_new, power_set_B2{j}(k, :))) = [ ];
                        [binary_decomp, span_B1] = R_is_span_union(B1_new, B2_new, R);
                        if binary_decomp == 1
                            
                            % Step 12
                            [model_N1, delta1] = deficiency(model, span_B1);
                            
                            % Step 13
                            if delta1 == 0
                                
                                % Step 14
                                PL_RDK = is_PL_RDK(model_N1, m);
                                if PL_RDK
                                    weakly_reversible = is_weakly_reversible(model_N1, m);
                                    if weakly_reversible
                                        reactant_complex1 = reactant_complex_list{SF_pair1};
                                        reactant_complex2 = reactant_complex_list{SF_pair2};
                                        same_linkage_class = in_same_linkage_class(reactant_complex1, reactant_complex2, model_N1, m);
                                        if same_linkage_class
                                            ACR_species{end+1} = [model.species{SF_pair(i, 3)} ' [SF-pair (R' num2str(SF_pair(i, 1)) ', R' num2str(SF_pair(i, 2)) ') | deficiency ' num2str(delta1) ']'];
                                            flag = 1;
                                            break
                                        end
                                    end
                                end
                            elseif delta1 == 1
                                
                                % Step 15
                                PL_RDK = is_PL_RDK(model_N1, m);
                                if PL_RDK
                                    reactant_complex1 = reactant_complex_list{SF_pair1};
                                    reactant_complex2 = reactant_complex_list{SF_pair2};
                                    if (is_nonterminal(reactant_complex1, model_N1, m) && is_nonterminal(reactant_complex2, model_N1, m))
                                        ACR_species{end+1} = [model.species{SF_pair(i, 3)} ' [SF-pair (R' num2str(SF_pair(i, 1)) ', R' num2str(SF_pair(i, 2)) ') | deficiency ' num2str(delta1) ']'];
                                        flag = 1;
                                        break
                                    end
                                end
                            end
                        end
                    end
                    if flag == 1
                        flag = 0;
                        break
                    end
                end
            end
        end
    end
    
    % Arrange the species in alphabetical order
    ACR_species = sort(ACR_species);
    
    
    
    %
    % Step 16: Display the results
    %
    
    % Case 1: No ACR in any species
    if numel(ACR_species) == 0
        disp(['The algorithm was not able to identify absolute concentration robustness in any species for ' model.id '.']);
        fprintf('\n\n')
    
    % Case 2: An ACR on a species exists
    else
        fprintf([model.id ' has absolute concentration robustness in the following species: \n\n']);
        fprintf('%s \n', ACR_species{:});
        fprintf('\n');
    end

end










% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% %                                                   % %
% % The following are functions used in the algorithm % %
% %                                                   % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                     %
% Function 1 of 11: model_species                                     %
%                                                                     %
%    - Purpose: To fill out list of species based on given reactions  %
%    - Input                                                          %
%         - model: empty species list                                 %
%    - Outputs                                                        %
%         - model: completed structure                                %
%         - m: number of species                                      %
%    - Used in                                                        %
%         - acr (Step 1)                                              %
%         - deficiency                                                %
%                                                                     %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [model, m] = model_species(model)

    % Initialize list of species
    model.species = { };

    % Get all species from reactants
    for i = 1:numel(model.reaction)
        for j = 1:numel(model.reaction(i).reactant)
            model.species{end+1} = model.reaction(i).reactant(j).species;
        end
    end
    
    % Get species from products
    for i = 1:numel(model.reaction)
        for j = 1:numel(model.reaction(i).product)
            model.species{end+1} = model.reaction(i).product(j).species;
        end
    end
    
    % Get only unique species
    model.species = unique(model.species);
    
    % Count the number of species
    m = numel(model.species);
    
end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                 %
% Function 2 of 11: stoich_matrix                                 %
%                                                                 %
%    - Purpose: To form the stoichometrix matrix N and the set of %
%          reactant and product complexes                         %
%    - Inputs                                                     %
%         - model: complete structure                             %
%         - m: number of species                                  %
%    - Outputs                                                    %
%         - N: stoichiometric matrix                              %
%         - reactant_complex: matrix of reactant complexes        %
%         - product_complex: matrix of product complexes          %
%         - r: total number of reactions                          %
%    - Used in                                                    %
%         - acr (Step 2)                                          %
%         - deficiency                                            %
%         - is_PL_RDK                                             %
%         - is_weakly_reversible                                  %
%         - in_same_linkage_class                                 %
%         - is_nonterminal                                        %
%                                                                 %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [N, reactant_complex, product_complex, r] = stoich_matrix(model, m)

    % Initialize the matrix of reactant complexes
    reactant_complex = [ ];
    
    % Initialize the matrix of product complexes
    product_complex = [ ];
    
    % Initialize the stoichiometric matrix
    N = [ ];
    
    % For each reaction in the model
    for i = 1:numel(model.reaction)
      
        % Initialize the vector for the reaction's reactant complex
        reactant_complex(:, end+1) = zeros(m, 1);
        
        % Fill it out with the stoichiometric coefficients of the species in the reactant complex
        for j = 1:numel(model.reaction(i).reactant)
            reactant_complex(find(strcmp(model.reaction(i).reactant(j).species, model.species), 1), end) = model.reaction(i).reactant(j).stoichiometry;
        end
        
        % Initialize the vector for the reaction's product complex
        product_complex(:, end+1) = zeros(m, 1);
        
        % Fill it out with the stoichiometric coefficients of the species in the product complex
        for j = 1:numel(model.reaction(i).product)
            product_complex(find(strcmp(model.reaction(i).product(j).species, model.species), 1), end) = model.reaction(i).product(j).stoichiometry;
        end
        
        % Create a vector for the stoichiometric matrix: Difference between the two previous vectors
        N(:, end+1) = product_complex(:, end) - reactant_complex(:, end);
        
        % If the reaction is reversible
        if model.reaction(i).reversible
          
            % Insert a new vector for the reactant complex: make it same as the product complex
            reactant_complex(:, end+1) = product_complex(:, end);
            
            % Insert a new vector for the product complex: make it the same as the reactant complex
            product_complex(:, end+1) = reactant_complex(:, end-1);
            
            % Insert a new vector in the stoichiometric matrix: make it the additive inverse of the vector formed earlier
            N(:, end+1) = -N(:, end);
        end
    end
    
    % Count the total number of reactions
    r = size(N, 2);

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                   %
% Function 3 of 11: kin_ord_matrix                                  %
%                                                                   %
%    - Purpose: To form the kinetic order matrix of the mass action %
%         system                                                    %
%    - Inputs                                                       %
%         - model: complete structure                               %
%         - m: number of species                                    %
%    - Output                                                       %
%         - F: kinetic order matrix                                 %
%    - Used in                                                      %
%         - acr (Step 3)                                            %
%         - is_PL_RDK                                               %
%                                                                   %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function F = kin_ord_matrix(model, m)
    
    % Initialize matrix F
    F = [ ];

    % Go through each reaction
    for i = 1:numel(model.reaction)
        
        % Case 1: The reaction is NOT reversible
        if model.reaction(i).reversible == 0
            
            % Add a row of zeros
            F(end+1, :) = zeros(1, m);
            
            % Fill out the kinetic order of all the species in the reactant
            for j = 1:numel(model.reaction(i).kinetic.reactant1)
                F(end, find(strcmp(model.reaction(i).reactant(j).species, model.species))) = model.reaction(i).kinetic.reactant1(j);
            end
        
        % Case 2: The reaction is reversible
        else
            
            % Add a row of zeros
            F(end+1, :) = zeros(1, m);
            
            % Fill out the kinetic order of all the species in the reactant in the first direction
            for j = 1:numel(model.reaction(i).kinetic.reactant1)
                F(end, find(strcmp(model.reaction(i).reactant(j).species, model.species))) = model.reaction(i).kinetic.reactant1(j);
            end
            
            % Add a row of zeros
            F(end+1, :) = zeros(1, m);
            
            % Fill out the kinetic order of all the species in the reactant in the other direction
            for j = 1:numel(model.reaction(i).kinetic.reactant2)
                F(end, find(strcmp(model.reaction(i).product(j).species, model.species))) = model.reaction(i).kinetic.reactant2(j);
            end
        end
    end

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                       %
% Function 4 of 11: R_is_span_union                                     %
%                                                                       %
%    - Purpose: To check if R is the union of the span(B1) and span(B2) %
%    - Inputs                                                           %
%         - B1: set of vectors to be extended to a basis                %
%         - B2: set of added vectors to form a basis                    %
%         - R: reaction matrix                                          %
%    - Ouputs                                                           %
%         - binary_decomp: logical; whether R is the union of span(B1)  %
%              and span(B2) or not                                      %
%         - span_B1: span(B1)                                           %
%    - Used in acr (STEPS 10, 11)                                       %
%                                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [binary_decomp, span_B1] = R_is_span_union(B1, B2, R)
    
    % Get the reaction vectors forming B1 and forming B2
    set_B1 = R(B1, :);
    set_B2 = R(B2, :);
    
    % Initialize list of reactions in span(B1) and in span(B2)
    span_B1 = [ ];
    span_B2 = [ ];
    
    % Go through each reaction vector
    for i = 1:size(R, 1)
        
        % If the reaction vector is a linear combination of set_B1
        if rank([set_B1', R(i, :)']) == rank(set_B1')
            
            % Add this reaction vector number to span_B1
            span_B1(end+1) = i;
        end
        
        % If the reaction vector is a linear combination of set_B2
        if rank([set_B2', R(i, :)']) == rank(set_B2')
            
            % Add this reaction vector number to span_B2
            span_B2(end+1) = i;
        end
    end
    
    % Get the union of span(B1) and span(B2)
    span_B1_U_span_B2 = union(span_B1, span_B2);
    
    % If the union becomes a vertical vector: This happens in Octave when one of the sets being combined is empty
    if size(span_B1_U_span_B2, 1) > 1
        
        % Transpose
        span_B1_U_span_B2 = span_B1_U_span_B2';
    end
    
    % Check if span_B1_U_span_B2 contains all reactions
    % If R is the union of span(B1) and span(B2), an independent binary decomposition is formed
    if isequal(1:size(R, 1), span_B1_U_span_B2)
        binary_decomp = 1;
    else
        binary_decomp = 0;
    end

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                       %
% Function 5 of 11: deficiency                                          %
%                                                                       %
%    - Purpose: To get the deficiency of a system created from span(B1) %
%    - Inputs                                                           %
%         - model: complete structure                                   %
%         - span_B1: span(B1)                                           %
%    - Ouputs                                                           %
%         - model_N1: new model created from span(B1)                   %
%         - delta1: deficiency of model_N1                              %
%    - Used in acr (Step 12)                                            %
%                                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [model_N1, delta1] = deficiency(model, span_B1)
    
    % Name the new model
    model_N1.id = 'span_B1';
    
    % Prepare the list of species
    model_N1.species = { };
    
    % Create a vector of model.reaction numbers for the total number of reactions
    reac_num = [ ];
    for i = 1:numel(model.reaction)
        if model.reaction(i).reversible == 0
            reac_num(end+1) = i;
        else
            reac_num(end+1) = i;
            reac_num(end+1) = i;
        end
    end
    
    % Get all model reaction numbers corresponding to the reactions in span_B1
    reac = unique(reac_num(span_B1));
    
    % Put these reaction numbers in model_N1
    for i = 1:numel(reac)
        model_N1.reaction(i) = model.reaction(reac(i));
    end
    
    % Add to model_N1.species all species indicated in the reactions of model_N1
    [model_N1, m] = model_species(model_N1);
    
    % Form stoichiometric matrix N
    [N, reactant_complex, product_complex, r] = stoich_matrix(model_N1, m);
        
    % Get just the unique complexes
    % index(i) is the index in all_complex of the reactant complex in reaction i
    [all_complex, ~, index] = unique([reactant_complex, product_complex]', 'rows');
    all_complex = all_complex';
    
    % Count the number of complexes
    n = size(all_complex, 2);
    
    % Initialize a matrix (complexes x complexes) for the reacts_to relation
    % This is for testing reversibility of the network
    reacts_to = false(n, n);
    
    % Initialize matrix (complexes x total reactions) for the reacts_in relation
    % This is the incidence matrix I_a
    reacts_in = zeros(n, r);
    
    % Fill out the entries of the matrices
    for i = 1:r
        
        % reacts_to(i, j) = true iff there is a reaction r: y_i -> y_j
        reacts_to(index(i), index(i + r)) = true;
        
        % reacts_in(i, r) = -1 and reacts_in(j, r) = 1) iff there is a reaction r: y_i -> y_j
        reacts_in(index(i), i) = -1;
        reacts_in(index(i+r), i) = 1;
    end
    
    % Linkage classes
    % Count number of connected components of an undirected graph
    linkage_class = conncomp(graph(reacts_to | reacts_to'));
    
    % Count the number of linkage classes
    l = max(linkage_class);
    
    % Get the rank of the reaction network
    % S = Im N
    % dim S = dim (Im N) = rank(N)
    % Note: We talk of "dimension of a linear transformation" and "rank of a matrix"
    s = rank(N);
    
    % Compute the deficiency of the reaction network
    delta1 = n - l - s;

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                 %
% Function 6 of 11: is_PL_RDK                                     %
%                                                                 %
%    - Purpose: To check if a system is a power law system with   %
%         reactant-determined kinetics (PL-RDK)                   %
%    - Inputs                                                     %
%         - model: complete structure                             %
%         - m: number of species                                  %
%    - Ouput                                                      %
%         - PL_RDK: logical; whether the system is PL-RDK or not  %
%    - Used in acr (STEPS 14, 15)                                 %
%                                                                 %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function PL_RDK = is_PL_RDK(model, m)
    
    % Get reactant_complex from function for stoichiometric matrix
    [~, reactant_complex, ~, ~] = stoich_matrix(model, m);
    
    % Form matrix of reactant complexes
    reactant_complex = reactant_complex';
    
    % Get the unique reactant complexes
    [reactant_complex_unique, ~, label] = unique(reactant_complex, 'rows');
    
    % Count the number of unique reactant complexes
    n_r = size(reactant_complex_unique, 1);
        
    % Get the labels of non-unique reactant complexes
    same_label = find(hist(label, unique(label)) > 1);
    
    % Initialize list reactions numbers with similar reactants
    branching_complex = { };
    
    % Group together reaction numbers with the same label
    for i = 1:numel(same_label)
        branching_complex{i} = find((label == same_label(i)));
    end
    
    % Initialize list of branching reactions
    branching_reaction = [ ];
    
    % Create a list of pairwise branching reactions
    for i = 1:numel(branching_complex)
        for j = 1:numel(branching_complex{i})
            for k = j+1:numel(branching_complex{i})
                
                % The list will have unique pairings: if (1, 2) is recorded, (2, 1) will no longer be noted
                branching_reaction(end+1, :) = [branching_complex{i}(j), branching_complex{i}(k)];
            end
        end
    end

    % Form kinetic order matrix F
    F = kin_ord_matrix(model, m);
    
    % Initialize list to track each pair
    F_identical = [ ];
    
    % Go through each pair
    for i = 1:size(branching_reaction, 1)
        
        % 1 means the pair has the same kinetic orders, 0 otherwise
        F_identical(i) = isequal(F(branching_reaction(i, 1), :), F(branching_reaction(i, 2), :));
    end
    
    % Case 1: There is at least 1 pair with different kinetic orders
    if ismember(0, F_identical)
        
        % The PLK system does NOT have reactant-determined kinetics (RDK)
        PL_RDK = 0;
    
    % Case 2: Each pair have the same kinetic order
    else
        
        % The PLK system has RDK
        PL_RDK = 1;
    end

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                     %
% Function 7 of 11: is_weakly_reversible                              %
%                                                                     %
%    - Purpose: To check if a system is weakly reversible             %
%    - Inputs                                                         %
%         - model: complete structure                                 %
%         - m: number of species                                      %
%    - Ouput                                                          %
%         - weakly_reversible: logical; whether the system is weakly  %
%              reversible or not                                      %
%    - Used in acr (Step 14)                                          %
%                                                                     %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function weakly_reversible = is_weakly_reversible(model, m)
    
    % Form stoichiometric matrix N
    [~, reactant_complex, product_complex, r] = stoich_matrix(model, m);
    
    % Get just the unique complexes
    % index(i) is the index in all_complex of the reactant complex in reaction i
    [all_complex, ~, index] = unique([reactant_complex, product_complex]', 'rows');
    all_complex = all_complex';
    
    % Count the number of complexes
    n = size(all_complex, 2);
        
    % Initialize a matrix (complexes x complexes) for the reacts_to relation
    % This is for testing reversibility of the network
    reacts_to = false(n, n);
    
    % Initialize matrix (complexes x total reactions) for the reacts_in relation
    % This is the incidence matrix I_a
    reacts_in = zeros(n, r);
    
    % Fill out the entries of the matrices
    for i = 1:r
        
        % reacts_to(i, j) = true iff there is a reaction r: y_i -> y_j
        reacts_to(index(i), index(i + r)) = true;
        
        % reacts_in(i, r) = -1 and reacts_in(j, r) = 1) iff there is a reaction r: y_i -> y_j
        reacts_in(index(i), i) = -1;
        reacts_in(index(i+r), i) = 1;
    end
    
    % Linkage classes
    % Count number of connected components of an undirected graph
    linkage_class = conncomp(graph(reacts_to | reacts_to'));
    
    % Count the number of linkage classes
    l = max(linkage_class);
        
    % Check if the network is reversibile
    is_reversible = isequal(reacts_to, reacts_to');
    
    % Strong linkage classes
    if is_reversible
        strong_linkage_class = linkage_class;
    else
        % Count number of connected components of a directed graph
        strong_linkage_class = conncomp(digraph(reacts_to));
    end
    
    % Count the number of strong linkage classes
    sl = max(strong_linkage_class);

    % Weakly reversible if the number of linkage classes and the number of strong linkage classes are the same
    if sl == l
        weakly_reversible = 1;
    else
        weakly_reversible = 0;
    end

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                               %
% Function 8 of 11: add_vertices                                %
%                                                               %
%    - Purpose: To add in a (di)graph all complexes as vertices %
%    - Inputs                                                   %
%         - g: (di)graph without vertices                       %
%         - n: number of complexes                              %
%         - Y: matrix of complexes                              %
%         - model: complete structure                           %
%    - Ouputs                                                   %
%         - g: (di)graph with all vertices                      %
%         - vertices: list of vertices in g                     %
%    - Used in                                                  %
%         - in_same_linkage_class                               %
%         - is_nonterminal                                      %
%                                                               %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [g, vertices] = add_vertices(g, n, all_complex, model)
    
    % Go through each column of Y (a complex)
    for i = 1:n
        
        % Get each column of Y
        all_complex_ = all_complex(:, i);

        % Get the index numbers of nonzero entries
        all_complex_nnz = find(all_complex_);

        % For the zero complex
        if numel(all_complex_nnz) == 0
            complex = '0';
        
        % Otherwise
        else
            
            % Check which species appear in the complex
            for j = 1:numel(all_complex_nnz)
                
                % For the first species
                if j == 1
                    
                    % Don't show the stoichiometry if it's 1
                    if all_complex_(all_complex_nnz(j)) == 1
                        complex = [model.species{all_complex_nnz(j)}];
                    
                    % Otherwise, include it
                    else
                        complex = [num2str(all_complex_(all_complex_nnz(j))), model.species{all_complex_nnz(j)}];
                    end
                
                % We need the + sign for succeeding species
                else
                    if all_complex_(all_complex_nnz(j)) == 1
                        complex = [complex, '+', model.species{all_complex_nnz(j)}];
                    else
                        complex = [complex, '+', num2str(all_complex_(all_complex_nnz(j))), model.species{all_complex_nnz(j)}];
                    end
                end
            end 
        end
        
        % Add this complex in the list of vertices of g
        g = addnode(g, complex);
    end
    
    % Convert table of nodes to a cell array
    vertices = table2cell(g.Nodes);

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                           %
% Function 9 of 11: add_edges                               %
%                                                           %
%    - Purpose: To add in a (di)graph all edges             %
%    - Inputs                                               %
%         - g: (di)graph without vertices                   %
%         - r: number of reactions                          %
%         - Y: matrix of complexes                          %
%         - reactant_complex: matrix of reactant complexes  %
%         - product_complex: matrix of product complexes    %
%    - Ouputs                                               %
%         - g: (di)graph with all vertices                  %
%         - edges: list of edges in g                       %
%    - Used in                                              %
%         - in_same_linkage_class                           %
%         - is_nonterminal                                  %
%                                                           %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [g, edges] = add_edges(g, vertices, r, Y, reactant_complex, product_complex)

    % Go through each reaction
    for i = 1:r

        % ~ suppresses the original output
        [~, loc1] = ismember(reactant_complex(:, i)', Y', 'rows');
        [~, loc2] = ismember(product_complex(:, i)', Y', 'rows');
        
        % Add edges to g: Ci -> Cj forms an edge
        g = addedge(g, vertices{loc1}, vertices{loc2});
    end

    % Convert table of edges to a cell array
    edges = table2cell(g.Edges);

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                         %
% Function 10 of 11: in_same_linkage_class                                %
%                                                                         %
%    - Purpose: To check if two complexes are in the same linkage class   %
%    - Inputs                                                             %
%         - complex1: first complex                                       %
%         - complex2: second complex                                      %
%         - model: complete structure                                     %
%         - m: number of species                                          %
%    - Ouput                                                              %
%         - same_linkage_class: logical; whether the two complexes are in %
%              the same linkage class or not                              %
%    - Used in acr (Step 14)                                              %
%                                                                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function same_linkage_class = in_same_linkage_class(complex1, complex2, model, m)
    
    % Form stoichiometric matrix N
    [~, reactant_complex, product_complex, r] = stoich_matrix(model, m);
    
    
    % Get just the unique complexes
    all_complex = unique([reactant_complex, product_complex]', 'rows');
    all_complex = all_complex';
    
    % Count the number of complexes
    n = size(all_complex, 2);
        
    % Initialize an undirected graph g
    g = graph();

    % Add to g all complexes as vertices
    [g, vertices] = add_vertices(g, n, all_complex, model);

    % Check if complex1 is a valid string of complexes
    if isempty(find(strcmp(vertices, complex1)))
        disp([complex1 ' is not a complex in the network.']);
    end
    
    % Check if complex2 is a valid string of complexes
    if isempty(find(strcmp(vertices, complex2)))
        disp([complex2 ' is not a complex in the network.']);
    end
    
    % If complex1 or complex2 is not valid, exit the function
    if (isempty(find(strcmp(vertices, complex1))) || isempty(find(strcmp(vertices, complex2))))
        same_linkage_class = [ ];
        return
    end
    
    % Add to g all reactions as edges
    g = add_edges(g, vertices, r, all_complex, reactant_complex, product_complex);
    
    % Get the linkage class of each complex    
    linkage_class = conncomp(g);
    
    % Determine if the given reactant complexes are in the same linkage class
    % 1 if they are in the same linkage class, 0 otherwise
    same_linkage_class = linkage_class(find(strcmp(vertices, complex1))) == linkage_class(find(strcmp(vertices, complex2)));

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                       %
% Function 11 of 11: is_nonterminal                                     %
%                                                                       %
%    - Purpose: To check if the complexes are nonterminal               %
%    - Inputs                                                           %
%         - check_complex: to check if nonterminal                      %
%         - model: complete structure                                   %
%         - m: number of species                                        %
%    - Ouput                                                            %
%         - nonterminal: logical; whether the complex is nonterminal or %
%              not                                                      %
%    - Used in acr (Step 15)                                            %
%                                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function nonterminal = is_nonterminal(check_complex, model, m)

    % Form stoichiometric matrix N
    [~, reactant_complex, product_complex, r] = stoich_matrix(model, m);
    
    % Get just the unique complexes
    all_complex = unique([reactant_complex, product_complex]', 'rows');
    all_complex = all_complex';
    
    % Count the number of complexes
    n = size(all_complex, 2);
        
    % Initialize a directed graph g
    g = digraph();
    
    % Add to g all complexes as vertices
    [g, vertices] = add_vertices(g, n, all_complex, model);
    
    % Check if check_complex is a valid string of complexes
    if isempty(find(strcmp(vertices, check_complex)))
        disp([check_complex ' is not a complex in the network.']);
        nonterminal = [ ];
        return
    end
    
    % Add a directed edge to g: Ci -> Cj forms an edge
    g = add_edges(g, vertices, r, all_complex, reactant_complex, product_complex);
    
    % Determine in which strong linkage class each complex belongs to
    strong_linkage_class = conncomp(g);
        
    % Initialize list of non-terminal strong linkage classes
    nonterminal_complexes = [ ];
    
    % Go through each complex
    for i = 1:numel(vertices)
        
        % Check to which complex it connects to (->)
        connections = successors(g, vertices{i});
        
        % If non-empty
        if isempty(connections) == 0
        
            % Go through each product complex
            for j = 1:numel(connections)
                
                % Check if the two complexes belong to the same strong linkage class
                same_strong_linkage_class = strong_linkage_class(i) == strong_linkage_class(find(strcmp(vertices, connections{j})));
                
                % If they do not, then the complex does not belong in a terminal strong linkage class
                if same_strong_linkage_class == 0
                    nonterminal_complexes(end+1) = i;
                end
            end
        end
    end
    
    % Generate the list of complexes that do not belong to a terminal strong linkage class
    nonterminal_complexes = unique(nonterminal_complexes);
    
    % Locate the other complexes belonging to the same strong linkage class as the ones identified to not belong in a terminal strong linkage class
    % These constitute the nonterminal complexes
    nonterminal_complexes = find(ismember(strong_linkage_class, strong_linkage_class(nonterminal_complexes)));
    
    % Check if the input complex is nonterminal
    if isempty(find(strcmp(vertices(nonterminal_complexes), check_complex))) == 1
        nonterminal = 0;
    else
        nonterminal = 1;
    end

end