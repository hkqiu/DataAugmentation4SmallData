There are the experimental codes for training, inference, visualization and data analysis.

# Part 1. Training
We have placed the training code in this file (train.ipynb) for the convenience of interested researchers to reproduce and use. 

Please note that we have tried multiple combinations of hyperparameters, but for the sake of brevity, we are showcasing the example with the best-performing set of hyperparameters.


# Part 2. Inference
We use the trained model to predict the glass transition temperature (Tg) of over 8,200,000 polyimide compounds. The prediction process can be found in the "inference.ipynb" file.

Please note that the prediction process may take a considerable amount of time, depending on the complexity of the model and the availability of computational resources. Therefore, before running "inference.ipynb," please ensure that your computing environment has sufficient computational power and time.

You can open "inference.ipynb" to view the prediction process.

# Part 3. Visualization
For detailed information and visualization regarding the various polymer datasets mentioned in **Figure 1b**, please refer to the "polymer_datasets_vis.ipynb" file.

The code for our proposed method, ACM (Atom Contribution Map), which visualizes fragment contribution rates and visualizes the model's inference mechanism, can be found in the "fragment_contribution.ipynb".

# Part 4. Data Analysis
We calculated the synthetic accessibility score (SAscore) for the candidate molecules to narrow down the screening range. The code for this part can be found in the "SAscore.ipynb" file.

In addition, we also calculated the number of rotatable bonds and the number of rings for the candidate molecules, which helps us gain deeper physical insights. The code for this part can be found in the "Rot2Ring.ipynb" file.


