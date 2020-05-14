# Model Predictive Control for Ackermann UGV
The scripts this repo entails are all the scripts and "scratch paper" related to my project in comparing different MPC algorithms for an ackermann-steered UGV drivepoint problem. It includes scripts for estimators of network parameters. All WIP.

## Overall Results
* A/B Beta Parameter estimation can be found in PlossEstimation.ipynb
* PID Control of time-delayed or high packet loss systems is not great
* Nonlinear prediction about as fast as Adaptive linear prediction, but Adaptive linear prediction significantly worse
* L-BFGS-B fastest optimizer for this, TNC/IPOPT most accurate, but SLSQP best performing. (Ipopt implemented with cyipopt, which can be found [here](https://github.com/matthias-k/cyipopt)
