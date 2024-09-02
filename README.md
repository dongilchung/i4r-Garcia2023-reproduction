## Instructions
The MATLAB code used to analyze the data presented in the article ```Reproduction report of Garcia et al. (2023) “Experiential values are underweighted in decisions involving symbolic options”``` can be found in this repository.

Each MATLAB script in the root folder begins by calling ```init.m```.

Before running any ```FigureX_TableX.m``` or ```Appendix_AX.m``` scripts, ensure that the ```data/fit``` folder contains the output from running ```fit_LE.m```, ```fit_logistic_ES.m```,  ```fit_logistic_EE.m```, ```fit_repl_LE.m```, ```fit_repl_logistic_ES.m```,  ```fit_repl_logistic_EE.m```, ```parameter_recovery_LE.m```, ```parameter_recovery_logistic_ES.m```, and  ```parameter_recovery_logistic_EE.m```.

The ```a_Run_Code.m``` script provides an overview of the order in which the code should be executed. Please refer to the explanation below for details on ```Hierarchical Bayesian Estimates```.

To get results of ```Hierarchical Bayesian Estimates```, ensure to run  ```trans_stanBehav_LE```, ```trans_stanBehav_ES```, and ```trans_stanBehav_EE``` and run R scripts in the ```Stan``` folder. 
(Scripts; 1. ```sv_mf_XX.R``` 2. ```Combine_compare_mf_XX.R```) 

For a successful computation, please run and combine at least 4 chains of each model.

## Data Structure
Behavioral data is located in the ```data/behavior/reformat/``` folder in CSV format.


| Variable    | Description                                                                                                     |
|-------------|-----------------------------------------------------------------------------------------------------------------|
| index       | Index identifier for each row                                                                                   |
| sub_id      | Subject identifier, identifying the participant                                                                 |
| phase       | Phase of the experiment (LE, ES, EE, SP)                                                                        |
| p1          | p(win) option 1                                                                                                 |
| p2          | p(win) option 2                                                                                                 |
| rtime       | Response time                                                                                                   |
| out         | Outcome of the choice (selected option)                                                                         |
| cfout       | Counterfactual outcome (unchosen option outcome)                                                               |
| cho         | Choice made by the participant (1 or 2)                                                                         |
| corr        | Correctness of the choice (1 if expected value of chosen option >= expected value of unchosen option, 0 else) |
| trial       | Trial number                                                                                                    |
| cond        | Condition of the experiment (3=60/40, 2=70/30, 1=80/20, 0=90/10)                                                 |
| chose_right | Whether the participant chose the rightmost option on screen                                                    |
| rew         | Total reward received (cumulated)                                                                               |
| sess        | Session identifier (When sess in (-1, -2), it means it was a training session, when sess is 0 or 1, it means first or second session) |
| op1         | Option 1 presented to the participant (filename or identifier)                                                   |
| op2         | Option 2 presented to the participant ((filename or identifier)                                                  |
| ev1         | Expected value of option 1                                                                                      |
| ev2         | Expected value of option 2                                                                                      |
| catch_trial | Indicates if it's a catch trial (testing attention)                                                              |
| reversed    | Indicates if the options were reversed in presentation (by default, option 1 is displayed on the left while option 2 is displayed on the right) |

