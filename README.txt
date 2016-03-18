                    #####################################
                    #    correlated_mutation_project    #
                    #####################################

                        ########################
                        #        Authors:      #
                        #     Andreu Bofill    #
                        #    Marina Reixachs   #
                        ########################


This python program studies correlated mutations in a single protein sequence or between two different protein sequences.

Given a single or two fasta files with the input sequences it runs Blast in order to find homologs to perform a MSA. Homologs are filtered according to Blast evalue in order that in the final MSA only one sequence for each species will be represented, the one with lower evalue. If two protein sequences are given, the species are matched in order to be contained the same species in both MSA.
With the MSA information a mutual information score is calculated as a measure of correlation between positions.
The scores for each position pair are provided in a tsv output file and also plotted as a heatmap that can be saved as a png file or visualised online in plotly.

Another input option of this program is to introduce a single or two multifasta files. If you want to avoid blast or you want to analyze an existing multifasta, you can. Both options fasta or multifasta are exclusive, you cannot use one fasta and one multifasta in the same execution. 

DEPENDENCIES:

The program runs with python3. Please check that it is installed:

	python3 --version

If it's not, you can install it the following way:

	sudo apt-get install python3

It also requires Biopython, Numpy and Matplotlib and Plotly.



PARAMETERS CONFIGURATION:

A parameters configuration file is provided. The main column of the file is invarable, and you need to add after tab space the needed information that in each line is given: The file has the following format:

	root	/root_of_the_future_output/
	blast	database	e-value_cutoff
	clustalw	/root_of_clustalw_program
	plotly	id	api_key
	Entrez_email	‘your_email_adress’


This option is optional, parameters.config file with a default options to be modified for each user and computer. You can also introduce a new file.config with the option (-c) but it needs to have the same format above.


USAGE:

usage: main.py [-h] [-i1 INFILE1] [-i2 INFILE2] [-mfa1 MULTIFASTA_1]
               [-mfa2 MULTIFASTA_2] [-c PARAMS] [-o OUTFILE] [-p HEATMAP]
               [-f FILTERED]

Correlated mutations

optional arguments:
  -h, --help            show this help message and exit
  -i1 INFILE1, --input1 INFILE1
                        Input file name
  -i2 INFILE2, --input2 INFILE2
                        Input file name
  -mfa1 MULTIFASTA_1, --multifasta1 MULTIFASTA_1
                        Multifasta Input file prepared to make a clustalW
                        alignment
  -mfa2 MULTIFASTA_2, --multifasta2 MULTIFASTA_2
                        Multifasta Input file prepared to make a clustalW
                        alignment
  -c PARAMS, --config PARAMS
                        Parameters configuration file
  -o OUTFILE, --output OUTFILE
                        Prefix for output files
  -p HEATMAP, --plot HEATMAP
                        Heatmap program options: png or plotly
  -f FILTERED, --filter FILTERED
                        If specified returns only positions with MI values
                        higher than cutt-off value


EXAMPLES:
	python3 main.py -i1 tfg-b.fa -o tfg-b_output

	python3 main.py -i1 tfg-b.fa -i2 tfg-b3.fa -p plotly

	python3 main.py -mfa1 multifasta1.mfa -p plotly -f 0.4

	python3 main.py -mfa1 multifasta1.mfa -mfa2 multifasta2.mfa -o output_name -f 0.04