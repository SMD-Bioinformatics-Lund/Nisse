#!/usr/bin/env python3

import sys
import gzip


def main():
    if len(sys.argv) != 3:
        print("Usage: <script> <in vcf (gz)> <cutoff>")
        sys.exit(1)
    
    in_fp = sys.argv[1]
    threshold = int(sys.argv[2])

    with gzip.open(in_fp, 'rt') as in_fh:
        for line in in_fh:
            line = line.rstrip()
            
            if line.startswith("#"):
                print(line)
                continue
    
            fields = line.split("\t")
            info_field = fields[7]
            key_vals = info_field.split(";")
            for key_val in key_vals:
                if "=" not in key_val:
                    continue
                (key, val) = key_val.split("=")
                if key == "RankScore":
                    score = float(val.split(":")[1])
       
     
        
                    if score >= threshold:
                        print(line)



if __name__ == "__main__":
    main()

