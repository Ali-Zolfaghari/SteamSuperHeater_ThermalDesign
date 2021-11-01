function [ data_flue ] = FlueProp( T,P,y_O2,y_N2,y_CO2,y_H2O )

Cp_O2  = refpropm('C','T',T+273.15,'P',100.0*P,'OXYGEN');
Rho_O2 = refpropm('D','T',T+273.15,'P',100.0*P,'OXYGEN');
Mu_O2  = refpropm('V','T',T+273.15,'P',100.0*P,'OXYGEN');
K_O2   = refpropm('L','T',T+273.15,'P',100.0*P,'OXYGEN');
Pr_O2  = refpropm('^','T',T+273.15,'P',100.0*P,'OXYGEN');

Cp_N2  = refpropm('C','T',T+273.15,'P',100.0*P,'NITROGEN');
Rho_N2 = refpropm('D','T',T+273.15,'P',100.0*P,'NITROGEN');
Mu_N2  = refpropm('V','T',T+273.15,'P',100.0*P,'NITROGEN');
K_N2   = refpropm('L','T',T+273.15,'P',100.0*P,'NITROGEN');
Pr_N2  = refpropm('^','T',T+273.15,'P',100.0*P,'NITROGEN');

Cp_CO2  = refpropm('C','T',T+273.15,'P',100.0*P,'CO2');
Rho_CO2 = refpropm('D','T',T+273.15,'P',100.0*P,'CO2');
Mu_CO2  = refpropm('V','T',T+273.15,'P',100.0*P,'CO2');
K_CO2   = refpropm('L','T',T+273.15,'P',100.0*P,'CO2');
Pr_CO2  = refpropm('^','T',T+273.15,'P',100.0*P,'CO2');

Cp_H2O  = refpropm('C','T',T+273.15,'P',100.0*P,'WATER');
Rho_H2O = refpropm('D','T',T+273.15,'P',100.0*P,'WATER');
Mu_H2O  = refpropm('V','T',T+273.15,'P',100.0*P,'WATER');
K_H2O   = refpropm('L','T',T+273.15,'P',100.0*P,'WATER');
Pr_H2O  = refpropm('^','T',T+273.15,'P',100.0*P,'WATER');

data_flue.cp  = y_O2*Cp_O2  + y_N2*Cp_N2  + y_CO2*Cp_CO2  + y_H2O*Cp_H2O ;
data_flue.rho = y_O2*Rho_O2 + y_N2*Rho_N2 + y_CO2*Rho_CO2 + y_H2O*Rho_H2O;
data_flue.mu  = y_O2*Mu_O2  + y_N2*Mu_N2  + y_CO2*Mu_CO2  + y_H2O*Mu_H2O ;
data_flue.k   = y_O2*K_O2   + y_N2*K_N2   + y_CO2*K_CO2   + y_H2O*K_H2O  ;
data_flue.pr  = y_O2*Pr_O2  + y_N2*Pr_N2  + y_CO2*Pr_CO2  + y_H2O*Pr_H2O ;

end
