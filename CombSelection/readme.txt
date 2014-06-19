*    selection.C: ROOT macro that actually performs selection

See file for further documentation.

*    conf.txt: configuration file that specifies sample name, cross section, flag for if is it in eos or not, and the base file path for each sample

For example, the first line is "HHToTTBB_14TeV 2.92 0 /afs/cern.ch/work/j/jlawhorn/public/HHToTTBB_14TeV/", which means it's going to look at files matching "/afs/cern.ch/work/j/jlawhorn/public/HHToTTBB_14TeV/HHToTTBB_14TeV*" using normal "ls" (the EOS/not flag is to toggle between "ls" and "eos ls", as well as toggle the "root://eoscms.cern.ch/" prefix.

You can comment out samples with a "#" at the beginning of the line.

*    master.sh: script that takes your conf.txt file and does pre-processing and optionally job submission

Run with "./master.sh conf.txt 1". For each line in conf.txt, generates a list of all the files in the sample (saved as <sample_name>.txt), computes the total number of events in the sample using getall.C (saved in <sample_name>_events.txt), then prepares a script to submit to Lxbtch (<sample_name>_run.sh) and submits the job (only if the last option is 1!)

getall.C: opens all ROOT files in a text file and finds the total number of events

