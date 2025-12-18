This repository contains R and Stan code for running the analysis published at https://dsts.dk/blog/2025-pl/. 

The model is a latent severity of signal model where each subject (or team in the application) has a latent severity (or strength in the application), which is estimated in a Bayesian framework. 

The model allows for the outcome to be time-to-event and random rater effects. 

This is a natural extension of the Bradley-Terry model type to handle time-to-event and random rater effects. The hierarchical Bayesian model structure handles these things quite easily.
