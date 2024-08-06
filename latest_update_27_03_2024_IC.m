%% Bcell Evolution 
function [NAb_ML]=latest_update_27_03_2024_IC(PI_Antibody_Aff,antigen_gen_var,ext_antibody_gen_var ,theta_c,delta_c)
%% Initialize all control variables
% Primordial controls


String_length=3;                                        % Sequence length
Num_alphabets=4;                                       % Number of alphabets = 4 in this case
NAb=1000;                                              % Initial number of B cells in the GC
Actual_NAb=NAb;                                        % Actual number of Bcells in the GC
Select_initial_MLforAb=1;                              % For initial conditions, select GC Bcells of a certain match length only
PI_Antibody_Af=PI_Antibody_Aff;
selffeedback=1;
% Germinal Center (GC) parameters, many from Luo and Perelson PNAS 2015
generation=1;                                          % Variable to keep track of B cell evolution. 
generation_max=100;                                    % Maximum number of generations
Ab_generation_time=1/2;                                % unit in days
Burstsize=4;                                           % Each B cell replicates twice
NAb_from_GC1=zeros(generation_max,1);                  % Number of B cells selected from a generation to burst, mutate and carry over to the next
mu_Ab=10^(-1);                                         % B cell mutation rate with a unit of /Bcell/generation,
                                                       % 10^-3/base/generation, order of 10^2 nucleotides
T_cellhelp=NAb/Burstsize;                              % Maximum number of cells that can be selected by T cells, Also prevents the GC from blowing up          
Avg_Ag=zeros(generation_max,1) ;                       % Average amount of antigen acquired by all the B cells in one generation
% Feedback parameters

    
% Array or matrix initializations for speedup
Bcell=uint8(zeros(NAb,String_length+2));            % Initialize Bcell array
Ab_selected=uint8(zeros(NAb,String_length+2));         % Initialize selected Bcell/antibodies array columns same as for Antibody matrix
ML=zeros(1,NAb);                                       % Initialize match lengths
Fit_T_Ab=zeros(1,NAb);                                 % Initialize T cell fitnesses
Virus=uint8(zeros(String_length,1));                   % Virus sequence length        
Avg_GC_Aff=zeros(generation_max,1);                     % Average match length of all the B cells in a single generation
Ag_BCR_I=100000;                                       % Total number of interactions of B cell with Ag
% Ag_BCR_IperBcell=Ag_interactions;                                  % Average Bcell - Ag interactions per B cell per generation 
prcnt_PC=0.05;                                         % prcnt of plasma cells
AvgPCML=zeros(generation_max,1);                       % Average feedback match length
ratio=zeros(generation_max,1);
N_max_PCR=NAb*prcnt_PC;                                % Maximum number of PCR cells in any given generation
PCR=zeros(generation,N_max_PCR)      ;                 % Plasma Cell repertoire
PCR_selectionindex=zeros(generation, N_max_PCR);       % Plasma Cell selection indices

NAb_ML_equals_3=zeros(generation_max,1);
NAb_ML_equals_2=zeros(generation_max,1);
NAb_ML_equals_1=zeros(generation_max,1);

IC_ML_1 = zeros(generation_max,1);
IC_ML_2 = zeros(generation_max,1);
IC_ML_3 = zeros(generation_max,1);


n_PC=zeros(generation_max,1);                          % matrix containing total number of plasma cells in each generation
n_MC=zeros(generation_max,1);                          % matrix containing total number of plasma cells in each generation

PC_ML_equals_3=zeros(generation_max,1);
PC_ML_equals_2=zeros(generation_max,1);
PC_ML_equals_1=zeros(generation_max,1);
Ag_Bcells=cell(generation_max,1) ;
Min_Ag=zeros(generation_max,1);
Max_Ag=zeros(generation_max,1);

antibody_cumm_ML_1 = zeros(generation_max,1);
antibody_cumm_ML_2 = zeros(generation_max,1);
antibody_cumm_ML_3 = zeros(generation_max,1);

PC_cumm_ML_1 = zeros(generation_max,1);
PC_cumm_ML_2 = zeros(generation_max,1);
PC_cumm_ML_3 = zeros(generation_max,1);

IC_cumm_aff_weighted = zeros(generation_max,1);

dummy_sum_ML_1 = 0;
dummy_sum_ML_2 = 0;
dummy_sum_ML_3 = 0;

dummy_sum_Ab_ML_1 = 0;
dummy_sum_Ab_ML_2 = 0;
dummy_sum_Ab_ML_3 = 0;

Ab_beta = 8.64*(10^7); % Parameters for endogenous antibodies
Ab_decay_rate = 0.01725;
PC_decay_rate = 0.015;

IC_cumm_aff_weighted_variable = 0;



eta_exo_max = 20; % scaling factor for antigen availability due to administered antibodies

 
Ag_steady = 8.5507*10^(4);
Ab_max = 0.06; % used in scaling to find antigen availability due to administered antibodies
antibody_excess_limit = 5; %antibody concentration beyond which they will not remain in excess.

conversion_mgkg_molar = 2.2222*(10^(-7));

avogadro_number = 6.022*(10^(23));

virus_number_molar = (1000/avogadro_number);
K_3bnc117 = (10^8);
k_ab1 = K_3bnc117/625; % rate constants for viral clearance by bNAb
k_ab2 = K_3bnc117/25;
k_ab3 = K_3bnc117;


IC_reference_exo = K_3bnc117*Ab_max*conversion_mgkg_molar*Ag_steady*virus_number_molar/(1 + K_3bnc117*Ab_max*conversion_mgkg_molar);

eta_reference_exo = 0.8*(eta_exo_max);
K_exo = IC_reference_exo*((eta_exo_max - eta_reference_exo)/eta_reference_exo);




% eta_reference_ini = 4;
eta_reference_ini = 2; % Antigen availability when no antibody therapy is done

IC_ini = ((eta_reference_ini/eta_exo_max)*K_exo)/(1-(eta_reference_ini/eta_exo_max));
Ab_ini = IC_ini/Ag_steady; % Antibody concentration when no antibody therapy is done
GCs_per_ml = 25; % Number of GCs per ml
avogadro_number = 6.022*(10^(23));


ratio_f = (eta_reference_ini/eta_exo_max)*(1/(k_ab1*Ag_steady*virus_number_molar));
Ab_low_before_pi = (ratio_f*K_exo)/(1 - ratio_f*K_exo*k_ab1 - ratio_f*k_ab1*Ag_steady*virus_number_molar); % Antibody concentration in molar when no antibody therpay is done



mol_wt_ab = 150000; % Molecular weight of bNAbs in dalton


dummy_conc_Ab_ML_1 = 0;
dummy_conc_Ab_ML_2 = 0;
dummy_conc_Ab_ML_3 = 0;
antigen_available = zeros(generation_max,1);
avg_FDC_aff = zeros(generation_max,1);



%% Initialization
 for k=1:String_length                                                           %% Generate random virus sequences
        Virus(k)=randi([1,Num_alphabets]);
 end
%% Antibody Feedback
            for j=1:NAb                                                          %% Generate initial Bcells/antibodies of a certain match length only
            
        while ( ML(j)~=Select_initial_MLforAb ) 
            for k=1:String_length
                Bcell(j,k)=randi([1,Num_alphabets]);
            end
            ML_forward=matchscore(Virus(:),Bcell(j,1:String_length));
            ML_reverse=matchscore(Virus(:),fliplr(Bcell(j,1:String_length)));
           ML(j)=max(ML_forward,ML_reverse);                                     %% Match length                                     
        end
            end
            
%% GC generation
for generation=1:generation_max
    
    exo_ab_ML_1 = 0;
    exo_ab_ML_2 = 0;
    exo_ab_ML_3 = 0;
    
    if PI_Antibody_Aff == 1
        exo_ab_ML_1 = ext_antibody_gen_var(generation);
    elseif PI_Antibody_Aff == 2
        exo_ab_ML_2 = ext_antibody_gen_var(generation);
    else
        exo_ab_ML_3 = ext_antibody_gen_var(generation);
    end
        
        
    total_ab_ML_1 = exo_ab_ML_1 + dummy_conc_Ab_ML_1;
    total_ab_ML_2 = exo_ab_ML_2 + dummy_conc_Ab_ML_2;
    total_ab_ML_3 = exo_ab_ML_3 + dummy_conc_Ab_ML_3;
    
    % Total Ab in molar concentration
    
    total_ab_molar_ML_1 = total_ab_ML_1*(conversion_mgkg_molar) + Ab_low_before_pi;
    total_ab_molar_ML_2 = total_ab_ML_2*(conversion_mgkg_molar);
    total_ab_molar_ML_3 = total_ab_ML_3*(conversion_mgkg_molar);
    
    
    antigen_conc(generation) = antigen_gen_var(generation)*virus_number_molar;
  
    IC_dinominator = 625 + K_3bnc117*(total_ab_molar_ML_1 + 25*total_ab_molar_ML_2 + 625*total_ab_molar_ML_3);
    
    

    IC_ML_1(generation) = K_3bnc117*antigen_conc(generation)*total_ab_molar_ML_1/IC_dinominator;
    IC_ML_2(generation) = 25*K_3bnc117*antigen_conc(generation)*total_ab_molar_ML_2/IC_dinominator;
    IC_ML_3(generation) = 625*K_3bnc117*antigen_conc(generation)*total_ab_molar_ML_3/IC_dinominator;
    
    IC_total(generation) = IC_ML_1(generation) + IC_ML_2(generation)+ IC_ML_3(generation); % IC concentrations for the particular generations
    
    
    
    
    
    
    Ag_BCR_IperBcell_total(generation) = eta_exo_max*(IC_total(generation)/(K_exo + IC_total(generation))); % Calculation of antigen availability using hill function

    antigen_available(generation) = Ag_BCR_IperBcell_total(generation);
    avg_FDC_aff(generation) = ((1/625)* IC_ML_1(generation) + (1/25)*IC_ML_2(generation) + IC_ML_3(generation))/(IC_ML_1(generation) + IC_ML_2(generation) + IC_ML_3(generation) );


        
        gen = generation;


            selection_index=zeros(1,NAb);                          % Initialize index of selected Abs, 3D array
            Fit_Ab=zeros(1,Actual_NAb);                                   % Initialize fitnesses
            array=zeros(1,Actual_NAb);                             % Dummy array used in selection of Bcells without replacement
            ratio(generation)=selffeedback;                                   % ratio of exogenous to endogenous ICs on FDC
%             
            for l=1:Actual_NAb       
            Bcell(j,String_length+1)=ML(j);    
            array(l)= l;                                           %Dummy matrix containing the position of B cell in a generation
            end               
          
            Avg_GC_Aff(generation)=mean(ML);                         
           

if (generation>2)

temp_PCR=zeros(1,n_PC(generation-2));                                
for i=1:n_PC(generation-2)                                          % LOop for extracting plasma cell matrix for the (current-2) generation (temp_PCR) from the net plasma cell matrix (PCR)
temp_PCR(i)=PCR((generation-2),i);
end
PCML=zeros(1,n_PC(generation-2));
PCML(:)=temp_PCR(:);                                                % PCML is dummy matrix for calculating average ML of plasma cells 
PCML(:)=[];                    
y=size(temp_PCR);                                                   % y is used to randomly pick plasma cells to assign antibodies to ICs
end
% LOop for calculating the amount of antigen acquired by all B cells at the end of AM
NO_SF=0;                                                            
NO_PI=0;                                                            
i=1;                                                                                                    


    while (i<=(Ag_BCR_IperBcell_total(generation)*Actual_NAb)) % Calculating proportion different affinity feedback antibodies on FDCs
        i = i+1;
        BCR=randi([1,Actual_NAb]);
        
        prob_IC_ML_1 = IC_ML_1(generation)/IC_total(generation);
        prob_IC_ML_2 = IC_ML_2(generation)/(IC_ML_2(generation) + IC_ML_3(generation));
        
        random_number_1 = rand;
        random_number_2 = rand;
        
       
        if random_number_1 < prob_IC_ML_1 
            IC_ML = 1;
            
        elseif random_number_2 < prob_IC_ML_2
            
            if (generation <3)
                IC_ML = 1;
            end
        
            if (generation >= 3)
                IC_ML = 2; 
            end
        else
            if (generation <3)
                IC_ML = 1;
            end
            
            if (generation >= 3)
                IC_ML = 3;
            end

        end
                    
                
                
                

            
            
        
        Fit_Ab(1,BCR)= (ML(1,BCR)-IC_ML +String_length)/(2*String_length); % Calculating probability of antigen acquisition
        temp_Fit_Ab=Fit_Ab(1,BCR);
        success1=0;
        r_MC=rand;
        if (r_MC<temp_Fit_Ab)
            success1=1;
        end
        if (success1==1) && (Bcell(BCR,String_length+2) < delta_c) 
            Bcell(BCR,String_length+2)=Bcell(BCR,String_length+2)+1;
        end
    end
Ag_Acquired=zeros(1,Actual_NAb);                                                    
Ag_Acquired(:)=Bcell(:,String_length+2);

Avg_Ag(generation)=mean(Ag_Acquired);   
Min_Ag(generation)=min(Ag_Acquired);                                 % Average amount of antigen acquired
Max_Ag(generation)=max(Ag_Acquired);                                 % Average amount of antigen acquired
Ag_Bcells(generation)={Ag_Acquired(:)}; 
j=0;        %  total number of B T cell interactions
k=1;        %  K compares with the population cap
             while((j<Actual_NAb) && (k<=T_cellhelp))                               	
               j=j+1;
                if isempty(array)==1
                break;
                end                                                                    % chose B cell at random without replacement
                dummy=size(array);    
                temp3=randi([1,dummy(2)]);
                temp2=array(temp3);                                                 % temp 3 stores the index of B cell on array and temp 2 is the identity of the B cell on the Antibody matrix
                
                if (double(min(Bcell(:,String_length+2))) == double(max(Bcell(:,String_length+2))))
                Fit_T_Ab(1,temp2) = 1;
                else
                Fit_T_Ab(1,temp2)= (double(Bcell(temp2,String_length+2))-double(min(Bcell(:,String_length+2))))/(double(max(Bcell(:,String_length+2)))-double(min(Bcell(:,String_length+2)))) ;  % Calculate fitnesss for B cell with the antibody complex
                end
                
                temp_Fit_Ab=Fit_T_Ab(temp2);
                success1=0;
                r_MC=rand;
                if (Bcell(temp2,String_length+2)<theta_c)
                elseif(r_MC<temp_Fit_Ab)
                    success1=1;
                end
                if (success1==1)                                                      % run selection MC move
                    selection_index(k)=temp2;                                          % Add B cell to list of selected cells with index stored in temp 2
                    k=k+1;
                    array(:,(temp3))=[];                                           % delete chpsen B cell index from array at temp 3 position
                end
                
             end  
             k=k-1;
             if (k==0||k==1)                                      % GC collapse termination loop 1
                for k=(generation+1):generation_max
                    NAb_ML_equals_3(k)= 0;
                    NAb_ML_equals_2(k)= 0;
                    NAb_ML_equals_1(k)= 0;                        
                end
                break;
              end
                n_PC(generation)=round((k-1)*prcnt_PC);                                 % number of plasma cells
                if(n_PC(generation)==0)
                n_PC(generation)=1;
                end
                for j=1:n_PC(generation)                                                 % Definining plasma cell repertoire for a generation
                X=randi([1,k]);
                PCR_selectionindex(generation,j)=selection_index(X);
                PCR(generation,j)=ML(selection_index(X));
                selection_index(X)=[];
                k=k-1;
                end

                for j=1:n_PC(generation)
                 
                if (PCR(generation,j)==3)                 
                   PC_ML_equals_3(generation) = PC_ML_equals_3(generation) + 1;
                elseif (PCR(generation,j)==2)                 
                   PC_ML_equals_2(generation) = PC_ML_equals_2(generation) + 1;
                elseif (PCR(generation,j)==1)
                   PC_ML_equals_1(generation) = PC_ML_equals_1(generation) + 1;
                end
                end
                
                dummy_sum_ML_1 = dummy_sum_ML_1*exp((-1)*PC_decay_rate) + PC_ML_equals_1(generation);
                dummy_sum_ML_2 = dummy_sum_ML_2*exp((-1)*PC_decay_rate) + PC_ML_equals_2(generation);
                dummy_sum_ML_3 = dummy_sum_ML_3*exp((-1)*PC_decay_rate) + PC_ML_equals_3(generation);
                
                
                PC_cumm_ML_1(generation) = dummy_sum_ML_1;
                PC_cumm_ML_2(generation) = dummy_sum_ML_2;
                PC_cumm_ML_3(generation) = dummy_sum_ML_3;
                
                dummy_sum_Ab_ML_1 = dummy_sum_Ab_ML_1*exp((-1)*Ab_decay_rate) + PC_cumm_ML_1(generation)*Ab_beta;
                dummy_sum_Ab_ML_2 = dummy_sum_Ab_ML_2*exp((-1)*Ab_decay_rate) + PC_cumm_ML_2(generation)*Ab_beta;
                dummy_sum_Ab_ML_3 = dummy_sum_Ab_ML_3*exp((-1)*Ab_decay_rate) + PC_cumm_ML_3(generation)*Ab_beta;
                
                dummy_conc_Ab_ML_1 = dummy_sum_Ab_ML_1*GCs_per_ml*(mol_wt_ab/avogadro_number)*(10^(3))*30; % Calculating cumulative antibody concentrations to calculating the feedback in the GCs 
                dummy_conc_Ab_ML_2 = dummy_sum_Ab_ML_2*GCs_per_ml*(mol_wt_ab/avogadro_number)*(10^(3))*30;
                dummy_conc_Ab_ML_3 = dummy_sum_Ab_ML_3*GCs_per_ml*(mol_wt_ab/avogadro_number)*(10^(3))*30;
                
                

                
                AvgPCML(generation)= (PC_ML_equals_3(generation)*3+PC_ML_equals_2(generation)*2+PC_ML_equals_1(generation))/(PC_ML_equals_3(generation)+PC_ML_equals_2(generation)+PC_ML_equals_1(generation));                                     % Average Plasma Cell Match Length is mean over the current generation
                n_MC(generation) = n_PC(generation);                                      % number of memory cells   
                if (n_PC(generation)==0)
                    n_MC(generation)=0;
                end
                for j=1:n_MC(generation)                                                  % Memory cell compartment selection and deletion 
                    X=randi([1,k]);
                    selection_index(X)=[];
                    k=k-1;
                end
                NAb_from_GC1(generation) = k;
                k=1;
               
                for h=1:NAb_from_GC1(generation)                                        % copying of selected B cells to dummy matrix
                    Ab_selected(k,:)=Bcell(selection_index(h),:);
                    k=k+1;
                end
                Bcell=uint8(zeros(NAb,String_length+2));                                 % reinitializing of B cell matrices
                ML=zeros(1,NAb);
                Actual_NAb=NAb_from_GC1(generation)*Burstsize;
                if (Actual_NAb==0 || Actual_NAb==1)                                      % GC collapse termination loop 2
                for k=(generation+1):generation_max                     
                      NAb_ML_equals_3(k)= 0;
                      NAb_ML_equals_2(k)= 0;
                      NAb_ML_equals_1(k)= 0;
                         
                end
                break;
                end
               k=1;
               for i=1:NAb_from_GC1(generation)                                         % Copying of B Cells to next generation
                    for j=1:Burstsize
                        Bcell(k,:)=Ab_selected(i,:);
                        k=k+1;
                    end
               end
               Bcell(k:NAb,:)=[];
               Bcell(1:(k-1),String_length+2)=0;
               ML(:,k:NAb)=[];
             
               count_mut=1;
              NAb_mu=round(mu_Ab*Actual_NAb);            % Number of mutated Abs per generation  
              while (count_mut~=(NAb_mu+1))       
                r_mut=randi([1,Actual_NAb]);     
              a=randi([1,String_length]);                    
              b=randi([1,Num_alphabets]);        
              if( Bcell(r_mut,a)~=b)
              Bcell(r_mut,a)=b;
              count_mut=count_mut+1;
              end
              end
            for j=1:Actual_NAb                    % determination of match lengths for new B cells which are to seed the next generation
                ML_forward=matchscore(Virus(:),Bcell(j,1:String_length));
                ML_reverse=matchscore(Virus(:),fliplr(Bcell(j,1:String_length)));
                ML(j)=max(ML_forward,ML_reverse); % Match length
                Bcell(j,String_length+1)=ML(j);
         
            if (ML(j)==3)                 
                   NAb_ML_equals_3(generation) = NAb_ML_equals_3(generation) + 1;
                elseif (ML(j)==2)                 
                   NAb_ML_equals_2(generation) = NAb_ML_equals_2(generation) + 1;
                elseif (ML(j)==1)
                   NAb_ML_equals_1(generation) = NAb_ML_equals_1(generation) + 1;
            end
            end

      
end     

        NAb_ML=zeros(generation_max,18); 

for i=1:generation_max
        NAb_ML(i,1)=NAb_ML_equals_1(i);
        NAb_ML(i,2)=NAb_ML_equals_2(i);
        NAb_ML(i,3)=NAb_ML_equals_3(i);
      
      
        NAb_ML(i,4)=NAb_from_GC1(i);
        NAb_ML(i,5)=Avg_GC_Aff(i);
        NAb_ML(i,6)=Avg_Ag(i);
        NAb_ML(i,7)=AvgPCML(i);
        NAb_ML(i,8)=ratio(i);

        NAb_ML(i,9)=PC_ML_equals_1(i);
        NAb_ML(i,10)=PC_ML_equals_2(i);
        NAb_ML(i,11)=PC_ML_equals_3(i);
       NAb_ML(i,12)=Min_Ag(i);
        NAb_ML(i,13)=Max_Ag(i);
        NAb_ML(i,14)=antigen_available(i);
        NAb_ML(i,15)=avg_FDC_aff(i);
        NAb_ML(i,16)= IC_ML_1(i);
        NAb_ML(i,17)= IC_ML_2(i);
        NAb_ML(i,18)= IC_ML_3(i);


       
end

end

 
function [ ms ] = matchscore(virus,antibody)
m=length(virus);
n=length(antibody);
ms=0;
L=zeros(m,n);

for i=1:m
    for j=1:n
        if virus(i)==antibody(j)
            if (i==1)||(j==1)
                L(i,j)=1;
            else
                L(i,j)=L(i-1,j-1)+1;
            end
            if L(i,j)>ms
                ms=L(i,j);
            end
        else
            L(i,j)=0;
        end
    end
end
end
