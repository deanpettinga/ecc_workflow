# Notes for Pierre's Pipeline
  1. save output of each step and ask program to look for partial output?
    - this could prevent having to re-run a long/costly step
    - start from partially processed reads?
    - if you successfully trimmed the reads, don't do it again.
  2. tell the program to create the temporary folders if they don't already exist.
    - my workflow would error out bc the tmp folders weren't recognized.
    - this means that i had to try to run the program, wait until it errored out bc the tmp path didn't exist, then i had to make the temp directory and rerun until i hit the same problem further down the line.
  3. 
