import shutil

def copy_template(yaml_template, yaml_new, job_template, job_new, temp, resname, num_iter, acc_name, time, jobname):
    """
    Creates a copy of a text file and performs a search-and-replace operation.

    :param file_path: Path to the original file.
    :param new_file_path: Path to the new file to be created.
    :param temp:
    :param resname:
    :param num_iter:
    """
    try:
        # Create a copy of the file
        shutil.copy(yaml_template, yaml_new)
        print(f"Copied '{yaml_template}' to '{yaml_new}'.")

        # Perform search and replace on the copied file
        with open(yaml_new, 'r') as file:
            content = file.read()

        content = content.replace("_TEMP_", temp)
        content = content.replace("_RESNAME_", resname)
        content = content.replace("_NUM_ITER_", num_iter)

        with open(yaml_new, 'w') as file:
            file.write(content)

        # Create a copy of the file
        shutil.copy(job_template, job_new)
        print(f"Copied '{job_template}' to '{job_new}'.")

        # Perform search and replace on the copied file
        with open(job_new, 'r') as file:
            content = file.read()

        content = content.replace("ACCOUNT_NAME", acc_name)
        content = content.replace("TIME", time)
        content = content.replace("JOBNAME", jobname)

        with open(job_new, 'w') as file:
            file.write(content)

    except Exception as e:
        print(f"An error occurred: {e}")

##  usage
yaml_template  = 'config_template.yaml'  # Replace with the path to your original text file
yaml_new = 'config.yaml'      # Replace with the desired path for the copied file
job_template = 'run_job_template.sh'
job_new = 'run_job.sh'

copy_template(yaml_template, yaml_new, job_template, job_new , temp="298.0", resname="MGLYOL", num_iter="10" , acc_name="energybio-eng", time="5:0:0", jobname="fep" )

