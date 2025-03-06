#!/usr/bin/env python3

import argparse
from pathlib import Path


DEST_HOST="rs-fs1.lunarc.lu.se"
PIPELINE_DEST="/fs1/pipelines/raredisease_rnaseq"
DEST="${DEST_HOST}:${PIPELINE_DEST}"

# echo "> Confirm, there are no running jobs for this pipeline? (y/n)"
# read -r no_running
# if ! [[ "${no_running}" =~ ^[yY]$ ]]; then
#     echo "Aborting"
#     exit 0
# fi

# Copy Nisse files
# bin files
# workflow files

# Copy Tomte files
# bin files
# workflow files

# Confirmations
# Are you on the master (main?) branch?


# Retrieve the current git.hash from each

# Copy config files

# Should this "deployer" actually live outside the pipelines?

def main(nisse_repo: Path, tomte_repo: Path, config_repo: Path):
    # Is master branch on latest commit?
    ...

def check_correct_branch():
    pass


def get_current_commit():
    pass


def parse_args():
    parser = argparse.ArgumentParser(description="Deploy Tomte and Nisse to the pipeline server")
    parser.add_argument("--nisse", help="Nisse directory", required=True)
    parser.add_argument("--tomte", help="Tomte directory", required=True)
    parser.add_argument("--configs", help="Config files repo", required=True)
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    main(args.nisse, args.tomte, args.configs)
