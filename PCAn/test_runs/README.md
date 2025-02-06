### Running PCAn on the test dataset

For testing PCAn, we provide _S. cerevisiae_ genome. After installing PCAn, please navigate to the test_runs directory to run this test.

1. Install PCAn as described [here](https://github.com/JHelsen/point-centromere-detection/tree/main/PCAn#readme)
2. Activate the conda PCAn environment and invoke PCAn
    ```bash
    conda activate pcan
    python /path/to/the/downloaded/PCAn/directory/AutomatedCENretrieval_ForPublishing.py
    ```
2. When you invoke PCAn, it will ask you for specific inputs (as outlined in the README for the PCAn folder). Please enter the inputs in the following order:
    ```   
     a. /path/to/the/downloaded/PCAn/directory/test_runs/Saccharomyces_cerevisiae.fna
    ```
    ```
     b. 0 (This selects the "Saccharomyces" genus to specify the centromere MEME motif to use for the run)
    ```
    ```
     c. Two options are available at this step:

       i. Yes (if you want to do the synteny checks using BLAST)

      ii. No  (if you do not want to do the synteny checks)
    ```

### Run Times
For this trial run, we ran PCAn on a Red Hat 8.5.0-22 Linux version 4.18.0-553.22.1.el8_10.x86_64 system. 

#### Without synteny checks
From the time we finished providing all the inputs, it took ~11 seconds to complete the run for the _S. cerevisiae_ genome
#### With synteny checks
From the time we finished providing all the inputs, it took ~121 seconds to complete the run for the _S. cerevisiae_ genome


### Sample outputs
The results of the PCAn run for the _S. cerevisiae_ genome are available in "sample_outputs". The BLAST folder is zipped to enable upload to the GitHub repository; it is normally a folder. More information on how to interpret the the PCAn outputs can be found [here](https://github.com/JHelsen/point-centromere-detection/tree/main/PCAn#pcan-outputs)
