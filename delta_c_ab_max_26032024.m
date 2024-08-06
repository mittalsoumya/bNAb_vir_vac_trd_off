clc;close all;clear all;
tic

for loop =  1:1:50
str = loop;
n = 1000; % Total days for which virus dynamics was run
gap_time = 0.001; % in days
time_span = 0:gap_time:n;
sz_time_span = size(time_span);
size_time_span = sz_time_span(2);


Max_PI_affinity = 3; % Maximum affinity of B cells
iterations = 100;
generation_max = 100; % One generation corresponds to 12 hours. One GC cycle is completeted in one generation.
Max_PI_Antibody_Aff = 3; % Maximum affinity of the administerd bNAbs
output_variables = 9 + 3*Max_PI_affinity;

theta_c = 2; % In the GC reaction
delta_c = 20;



beta = 2.4*10^(-8);  
delta = 1;
p = 978;
c = 23;
lambda_1 = 10^(5);
d_t = 0.1;
k_clearance = 0.0345; % Parameters for virus dynamics system of equation
Ag_steady = ((((beta*lambda_1*p)/(delta*d_t*c)) - 1)*(d_t/beta));
Infected_cell_steady = ((((beta*lambda_1*p)/(delta*d_t*c)) - 1)*((d_t*c)/(beta*p)));
Uninfected_cell_steady  = ((delta*c)/(beta*p));



mug_per_ml_mg_per_kg = (30/1000);
Antibody_input_1 = 0:0.1:10;  Antibody_input_2 = 15:5:1200; % Dose of administered bNAbs in (mug/ml)
Antibody_input = [Antibody_input_1 Antibody_input_2];


%% We calculate the bNAb on the viral load
Size_Antibody_input = size(Antibody_input);
Size_Ab_input = Size_Antibody_input(2);
GCrealization_iter = zeros(generation_max,output_variables, iterations,Max_PI_Antibody_Aff,Size_Ab_input);
for ab_input_index = 1:1:Size_Ab_input
    Ab_input = Antibody_input(ab_input_index);
    
    alpha_1 = 0.7809267*Ab_input;% Parameters for bNAbs
    alpha_2 = 0.2190733*Ab_input;
    eta_1 = 1.0452;
    eta_2 = 0.0945;
    
    for PI_Affinity = 1:1:Max_PI_Antibody_Aff
        
      
        k1 = (25^(PI_Affinity-1))*(k_clearance/(25^(Max_PI_Antibody_Aff-1)));
        
        y0 = [Ag_steady;Infected_cell_steady;Uninfected_cell_steady];
        F = @(t,y) [-k1*y(1)*(alpha_1*exp(-eta_1.*t) + alpha_2*exp(-eta_2.*t)) + p*y(2) - c*y(1); beta*y(1)*y(3) - delta*y(2); lambda_1 - d_t*y(3) - beta*y(1)*y(3) ];
        [t,y] = ode45(F, time_span, y0);

        
        for pp = 1:1:size_time_span

          antigen_var(pp) = y(pp,1); % Virus concentration with time
        end
        
        for pp = 1:1:size_time_span
            tt = (pp - 1)*gap_time;
            ext_antibody_var(pp) = (alpha_1*exp(-eta_1*tt) + alpha_2*exp(-eta_2*tt))*mug_per_ml_mg_per_kg; % bNAb decay with time
        end
        
        
        for i = 1:1:1001
            antigen_gen_var(i) = antigen_var(1/(2*gap_time)*(i-1)+1);
        end
        
        for i = 1:1:1001
            ext_antibody_gen_var(i) = ext_antibody_var(1/(2*gap_time)*(i-1)+1);% bANb concentration with time
        end
        


        antigen_gen_var = antigen_gen_var;
        ext_antibody_gen_var = ext_antibody_gen_var;
        
        
        
        parfor iter = 1:1:iterations
             
            GCrealization=latest_update_27_03_2024_IC(PI_Affinity,antigen_gen_var ,ext_antibody_gen_var, theta_c, delta_c );

            GCrealization_iter(:,:,iter,PI_Affinity,ab_input_index) = GCrealization(:,:); % Coupling between the virus dynamics and GC reaction
            
        end
    end
end


str  = num2str(loop);
string = ['vir_vac_hum_par_' str];
save(string)

end
toc
