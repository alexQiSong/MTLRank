import argparse
import os
import subprocess
from preprocess_expressions import get_expressions
from preprocess_get_tf_activities import get_tf_activities

def main():
    
    # Get argument
    parser = argparse.ArgumentParser(
                        prog = 'download_sample_data_and_preprocess',
                        description = 'This script will download sample data sets for liver and spleen tissue from HuBMAP data portal and preprocess the data and format it as input for MTLRank model training and GRN inference.',
                        epilog = '')
    parser.add_argument("-n","--n_jobs",default = 2)
    args = parser.parse_args()
    n_jobs = int(args.n_jobs)
    
    # Run preprocessing steps.
    get_expressions()
    get_tf_activities(n_jobs = n_jobs)

    return None

if __name__ == "__main__":
    main()