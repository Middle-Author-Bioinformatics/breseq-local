#!/bin/bash
eval "$(/home/ark/miniconda3/bin/conda shell.bash hook)"
conda activate base  # Activate the base environment where `boto3` is installed

exec > >(tee -i /home/ark/MAB/breseq/breseq_looper.log)
exec 2>&1

## Debugging information
#echo "Script started at $(date)"
#echo "Current directory: $(pwd)"
#echo "Environment variables: $(env)"

eval "$(/home/ark/miniconda3/bin/conda shell.bash hook)"
conda activate base  # Activate the base environment where `boto3` is installed

# Debug PATH
#echo "Updated PATH: $PATH"

# Debug Python environment
#echo "Python version being used:"
#which python3
#python3 --version
#
#echo "Python modules installed:"
#python3 -m pip list

KEY=$1
ID=$KEY
DIR=/home/ark/MAB/breseq/${ID}
OUT=/home/ark/MAB/breseq/completed/${ID}


name=$(grep 'Name' ${DIR}/form-data.txt | cut -d ' ' -f2)
email=$(grep 'Email' ${DIR}/form-data.txt | cut -d ' ' -f2)
platform=$(grep 'Platform' ${DIR}/form-data.txt | cut -d ' ' -f3)
ref=$(grep 'Reference' ${DIR}/form-data.txt | cut -d ' ' -f2)
mode=$(grep 'polymorphism-prediction' ${DIR}/form-data.txt | cut -d ' ' -f2)
metadata=$(grep 'Metadata' ${DIR}/form-data.txt | cut -d ' ' -f3)
# code to check if mode is empty
if [ -z "$mode" ]; then
    mode="clone"
fi

# Verify email
result=$(python3 /home/ark/MAB/bin/breseq-local/check_email.py --email ${email})
echo $result

# Set PATH to include Conda and script locations
export PATH="/home/ark/miniconda3/bin:/usr/local/bin:/usr/bin:/bin:/home/ark/MAB/bin/breseq-local:$PATH"
eval "$(/home/ark/miniconda3/bin/conda shell.bash hook)"
conda activate breseq_env

if [ $? -ne 0 ]; then
    echo "Error: Failed to activate Conda environment."
    exit 1
fi
sleep 5

# Run Breseq
# **************************************************************************************************
mkdir -p /home/ark/MAB/breseq/completed/${ID}
/home/ark/MAB/bin/breseq-local/generate_commands_for_samples.sh ${DIR}/${metadata} ${DIR}/form-data.txt ${DIR}/commands ${platform} ${ID}
for file in ${DIR}/commands/command-*; do
    bash $file
done
# **************************************************************************************************
mkdir -p ${OUT}/Combined_Summary
mkdir -p ${OUT}/analysis
for i in ${OUT}/breseq_*; do
    /home/ark/MAB/bin/breseq-local/breseq_parser.py -b ${i} -o ${OUT}/Combined_Summary -gbk ${DIR}/${ref} -m ${mode}
    /home/ark/MAB/bin/breseq-local/breseq_plotter.py ${i} ${OUT}/analysis
done
# **************************************************************************************************
/home/ark/MAB/bin/breseq-local/breseq_combiner.py -i ${OUT}/Combined_Summary -m ${mode}
# **************************************************************************************************
# **************************************************************************************************
# **************************************************************************************************
if [ $? -ne 0 ]; then
    echo "Error: Breseq failed."
    conda deactivate
    exit 1
fi
conda deactivate
sleep 5

# Removing misc output files
mv ${OUT} /home/ark/MAB/breseq/tmp/${ID}
mkdir /home/ark/MAB/breseq/completed/${ID}
for i in /home/ark/MAB/breseq/tmp/${ID}/breseq_*; do
    sample=$(basename $i)
    mv $i/output /home/ark/MAB/breseq/completed/${ID}/${sample}
done

mv /home/ark/MAB/breseq/tmp/${ID}/Combined_Summary /home/ark/MAB/breseq/completed/${ID}/
#rm -rf /home/ark/MAB/breseq/tmp/${ID}

# Archive results
tar -cf /home/ark/MAB/breseq/completed/${ID}.tar /home/ark/MAB/breseq/completed/${ID} && gzip /home/ark/MAB/breseq/completed/${ID}.tar

# Upload results to S3 and generate presigned URL
results_tar="/home/ark/MAB/breseq/completed/${ID}.tar.gz"
s3_key="${ID}.tar.gz"
python3 /home/ark/MAB/bin/breseq-local/push.py --bucket binfo-dump --output_key ${s3_key} --source ${results_tar}
url=$(python3 /home/ark/MAB/bin/breseq-local/gen_presign_url.py --bucket binfo-dump --key ${s3_key} --expiration 86400)

# Send email
python3 /home/ark/MAB/bin/breseq-local/send_email.py \
    --sender ark@midauthorbio.com \
    --recipient ${email} \
    --subject "Your Breseq Results!" \
    --body "Hi ${name},

    Your Breseq results are available for download using the link below. The link will expire in 24 hours.

    ${url}

    Please visit https://barricklab.org/twiki/pub/Lab/ToolsBacterialGenomeResequencing/documentation/output.html for documentation.

    Please reach out to ark@midauthorbio.com, or send us a note on https://midauthorbio.com/#contact if you have any questions.

    Thanks!
    MAB Team"

if [ $? -ne 0 ]; then
    echo "Error: send_email.py failed."
    conda deactivate
    exit 1
fi

sleep 5

#sudo rm -rf ${DIR}

conda deactivate
echo "Breseq completed successfully."



