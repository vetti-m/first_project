#!/bin/bash
echo "Running the script to get Metazoa Data..."
echo "Preparing files and folders..."

# run preparation bash script
bash prepare_pipeline.sh

echo "Successfully prepared, now starting the Metazoa Pipeline."
echo "Logging will be saved in doc/metazoa_pipeline.log."

# be sure conda is activated
eval "$(conda shell.bash hook)"
conda activate base

# run python script
python3 metazoa_pipeline.py


# check if the python script ran successfully
if [ $? -eq 0 ]; then        # check if the exit code is 0 (success)
    echo "-------------------------------------------------------------"
    echo "Metazoa Pipeline completed successfully :)"
    echo "Results are saved in the results directory."
    echo "-------------------------------------------------------------"
else
    echo "-------------------------------------------------------------"
    echo "Metazoa Pipeline encountered an error :("
    echo "Please check doc/metazoa_pipeline.log for more details."
    echo "-------------------------------------------------------------"
    exit 1     # exit with error code
fi