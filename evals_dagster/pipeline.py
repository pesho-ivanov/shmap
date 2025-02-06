from dagster import job, op, Field, In, Out
import subprocess
import os

@op(
    config_schema={
        "shmap_bin": Field(str, description="Path to the SHMAP binary."),
        "ref": Field(str, description="Path to the reference file."),
        "reads": Field(str, description="Path to the reads file."),
        "one_read": Field(str, description="Path to a single read file."),
        "output_prefix": Field(str, description="Prefix for output files."),
        "k": Field(int, description="Value of k parameter."),
        "r": Field(float, description="Value of r parameter."),
        "t": Field(float, description="Value of t parameter."),
        "min_diff": Field(float, description="Minimum difference parameter."),
        "max_overlap": Field(float, description="Maximum overlap parameter."),
        "metric": Field(str, description="Metric for evaluation."),
        "shmap_args": Field(str, default_value="", description="Additional SHMAP arguments.")
    },
    out=Out(str)
)
def run_shmap(context):
    config = context.op_config
    
    # Ensure output directory exists
    output_dir = os.path.dirname(config["output_prefix"])
    os.makedirs(output_dir, exist_ok=True)
    
    index_time_file = f"{config['output_prefix']}.index.time"
    time_file = f"{config['output_prefix']}.time"
    log_file = f"{config['output_prefix']}.log"
    paf_file = f"{config['output_prefix']}.paf"

    # Step 1: Indexing (if needed)
    context.log.info("Starting SHMAP indexing...")
    indexing_command = [
        config["shmap_bin"], "-s", config["ref"], "-p", config["one_read"],
        "-k", str(config["k"]), "-r", str(config["r"]),
        "-t", str(config["t"]), "-d", str(config["min_diff"]),
        "-o", str(config["max_overlap"]), "-m", config["metric"]
    ]
    if config["shmap_args"]:
        indexing_command.extend(config["shmap_args"].split())
    
    # Run indexing command and time it
    with open(index_time_file, 'w') as index_time_out:
        subprocess.run(["/usr/bin/time", "-f", "%U\t%M"] + indexing_command, stdout=subprocess.DEVNULL, stderr=index_time_out)
    
    # Step 2: Alignment
    context.log.info("Starting SHMAP alignment...")
    alignment_command = [
        config["shmap_bin"], "-s", config["ref"], "-p", config["reads"],
        "-z", f"{config['output_prefix']}.params", "-k", str(config["k"]),
        "-r", str(config["r"]), "-t", str(config["t"]),
        "-d", str(config["min_diff"]), "-o", str(config["max_overlap"]),
        "-m", config["metric"]
    ]
    if config["shmap_args"]:
        alignment_command.extend(config["shmap_args"].split())

    # Run alignment command and log output
    with open(time_file, 'w') as time_out, open(log_file, 'w') as log_out:
        subprocess.run(["/usr/bin/time", "-f", "%U\t%M"] + alignment_command, stdout=open(paf_file, 'w'), stderr=log_out)

    context.log.info("SHMAP evaluation completed.")
    return paf_file

@job
def shmap_job():
    run_shmap()
