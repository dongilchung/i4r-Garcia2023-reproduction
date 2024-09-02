clc
clear all
close all

%% 1) Author's code replication
fit_LE
fit_logistic_ES
fit_logistic_EE

%% 2) MLE replication
fit_repl_LE
fit_repl_logistic_ES
fit_repl_logistic_EE

%% 3) Parameter recovery
parameter_recovery_LE
parameter_recovery_logistic_ES
parameter_recovery_logistic_EE

%% 4) Transform mat data form to stan data form
trans_stanBehav_LE
trans_stanBehav_ES
trans_stanBehav_EE

%% 5) Figure-Table replication
Figure1_Table1
Figure2_Table2
Figure3_Table3
Figure4_Table4

%% 6) Appendix replication
Appendix_A1
Appendix_A2
Appendix_A3
Appendix_A4
Appendix_A5
Appendix_A6 % long time
Appendix_A7
Appendix_A8
Appendix_A9
Appendix_A10
Appendix_A11
Appendix_A12