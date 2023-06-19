[![Contributors][contributors-shield]][contributors-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]

# drugdiscoverychallenge
This resporitory contains new ideas about implementing Machine Learning for ALDH1 inhibition. ALDH1 is a well-known enzyme that contributes in various diseases, including cancer, Alzheimer's disease, and cardiovascular disease. This approach of finding new drugs is called Computer Asisted Drug Discovery (CADD). <br>
This specific resporitory uses a DecisionTree model to find 100 new Inhibitors for ALDH1. The final code that was obtained is visible in <b>main.ipynb</b>.

## DecisionTree model 
For the descision tree model SKlearn was used. First the 14 best descriptors were selecthed through a pre-computed PCA. Then a pipeline was initiated containing a MinMaxScaler, PCA and Decision tree classifier. No further settings were applied to this these components of the pipeline. The pipeline was fitted with a training set of about 1600 known inhibitors. 400 additional inhibitors were used to test the model. An example from the decision tree can be seen [here](output/descision_tree.pdf). 


## Tanimoto similarity
After the model was calibrated to reasonable scores, the model was used on a database set containing 10.000 Inhibitors from which the effect on ALDH1 is unknown. All Inhibitors which are working according to the model were then extracted. These Inhibitors were tested for Tanimoto similarity with the dataset from which the Inhibition effectiveness is known. Finally, 100 Inhibitors with the highest similarity were selected. <br>

## How the code works
This code runs by simply opening the <b>main.ipynb</b> notebook and running all lines. The code will produce a machine learning model and it will create some files which can be seen in the output folder. e.g. a decision_tree.pdf is made and some score figures are produced. 

<br>
Please see the report for more insights into this research. 

<br><i>For furhter questions please contact group 10</i>

[contributors-shield]: https://img.shields.io/github/contributors/Undefined1234/drugdiscoverychallenge.svg?style=for-the-badge
[contributors-url]: https://github.com/Undefined1234/drugdiscoverychallenge/contributors
[stars-shield]: https://img.shields.io/github/stars/Undefined1234/drugdiscoverychallenge.svg?style=for-the-badge
[stars-url]: https://github.com/Undefined1234/drugdiscoverychallenge/stargazers
[issues-shield]: https://img.shields.io/github/issues/Undefined1234/drugdiscoverychallenge.svg?style=for-the-badge
[issues-url]: https://github.com/Undefined1234/drugdiscoverychallenge/issues