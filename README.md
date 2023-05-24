## Full Sample Scripts for Running Simulation Study

These files are intended to run the Coordinated Multi-Neighborhood Learning (CML) and Single Neighborhood Learning (SNL) algorithms using the full sample version of these algorithms, including Markov Blanket estimation for a simulation study.

In order to use, first install the CML package and then download this repository to your machine.

Then, download the rds files for the networks you are interested in from the [bnlearn repository](https://www.bnlearn.com/bnrepository/).

Adjust the file names in "initialize.R" and "arrayscript.R" to match where you are storing this folder and where the results of the simulation study should be stored. Pay careful attention to where the RDS files should be placed for the networks.

Run the "initialize.R" to set up the results folder and the settings for each simulation (stored in "sim_vals.csv")

Then, you may run "arrayscript.R", which will run the setting defined by the first row in "sim_vals.csv". You can vary the array_number variable in "arrayscript.R" to change which setting will be run.
