# Notes for Pierre's Pipeline
  1. save output of each step and ask program to look for partial output?
    - this could prevent having to re-run a long/costly step
    - start from partially processed reads?
    - if you successfully trimmed the reads, don't do it again.
  2. tell the program to create the temporary folders if they don't already exist.
    - my workflow would error out bc the tmp folders weren't recognized.
    - this means that i had to try to run the program, wait until it errored out bc the tmp path didn't exist, then i had to make the temp directory and rerun until i hit the same problem further down the line.
  3. output should always have the sample as a prefix and the name of the file as the suffix.
    - i had an issue due naming and creating my files in different directories.
    - when i tried to go from `generate_bam_file.sh` to `call_ecc_regions`, ecc_caller could not find the appropriate bams because it assumed they were in the base directory.
    - `generate_bam_file` places -s before a file suffix: `{-s}.sorted.mergedandpe.bwamem.bam`
    - `call_ecc_regions` places -s in the middle of a file suffix `filtered.sorted.{-s}.bam`
  4. error messages:
    - is it possible to make error messages refer to a line in a script (bash or python?)
    - make the program stop if certain files are empty and give an error message. this would make my workflow error out, stop and give me a reason why.
  
