"""
Must Have ANVIO Installed
	Must be in ANVIO's Enviroment. If on Terminal use: $ conda activate anvio-6.2

Other Packages You Need To Install First:
	$ anvi-setup-ncbi-cogs     # For Genome Annotation in Make Storage
	$ conda install bowtie2
	$ conda install sqlite
	$ conda install biopython  # For Editing Fasta Files to Remove Non-ATGCN Character in Data Preparation
	$ conda install -c anaconda matplotlib ( conda install -c conda-forge/label/testing matplotlib )
	$ conda install illumina-utils

The following code is based on the ANVIO outlines:
	MetaPangenome: http://merenlab.org/data/prochlorococcus-metapangenome/#recruiting-and-profiling-reads-form-metagenomes
	Quality Filtering: http://merenlab.org/tutorials/assembly-based-metagenomics/#quality-filtering
	SNV: http://merenlab.org/2015/07/20/analyzing-variability/#generating-a-snv-scv-or-saav-profile
	THe SNV Part of: http://merenlab.org/data/sar11-saavs/

To Visualize SNVs (Not Used in This Code):
	Type R in the command line. Once there types the following: install.packages(c('ggplot2', 'reshape2', 'reshape', 'gridExtra', 'grid', 'plyr', 'gtools'))
"""

# --------------------------------------IMPORTS--------------------------------------------------------------------------#


# Modules Always Needed
import os
# Modules for ANVIO Storage
import multiprocessing
import csv
# Modules Needed to Generate Collection
import sqlite3
# To Edit FASTA Files and Remove Extended Nucleotides
from Bio.SeqRecord import SeqRecord
from Bio import Alphabet
from Bio.Seq import Seq
from Bio import SeqIO
# Get Inports from Command Line
import sys


# ----------------------------------------------------------------------------------------------------------------#

class helper_Files:
	def __init__(self, output_File_Directory):
		# Make Folder for Pangenome Files if not present. If it already exists, this wont do anything
		self.fasta_Files =  output_File_Directory + "fasta_Files/"
		self.contigs_Database = output_File_Directory + "contigs_Database/"
		self.filter_Directory = output_File_Directory + "Quality_Filtering_Files/"
		self.bam_Directory = output_File_Directory + "Bam_Files/"
		self.collection_Files =  output_File_Directory + "collection_Files/"
		self.output_SNV = output_File_Directory + "SNV_Files/"

	def make_Output_Folders(self, output_File_Directory):
		# Make Folder if not present. If it already exists, this wont do anything
		self.execute_Command("mkdir " + output_File_Directory)
		self.execute_Command("mkdir " + self.fasta_Files)
		self.execute_Command("mkdir " + self.contigs_Database)
		self.execute_Command("mkdir " + self.filter_Directory)
		self.execute_Command("mkdir " + self.bam_Directory)
		self.execute_Command("mkdir " + self.collection_Files)
		self.execute_Command("mkdir " + self.output_SNV)

	def execute_Command(self, command):
		print("\n Executing the Following Command: ", command)
		os.system(command)


class user_Inputs:
	def gather_Inputs(self):
		# Record Variables if Inputed Through Command Line. Example: 'python Run_Full_SNV.py Hello Sam' Translates to sys.argv = ['Run_Full_SNV.py', 'Hello', 'Sam']
		# Base Case: No Arguments
		if len(sys.argv) == 1:
			return "", "", "", "", "", "", "", "", "", "", ""  # Return Blanks (Using Names in 'if __name__ == "__main__"' at bottom
		# User Inputed Arguments
		print("Program Will Assume you Input Information as:")
		print("python Run_Full_Pangenome.py -in [Input File Directory] -out [Output File Directory] -S [Fasta File of Interest] -bin [Bin of Interest] -name [Project Name] -title [SNV Title] -T [Threads] -type [Profile Type: NT, AA, CDN] --full_run --profile_SCVs --low_storage \n")
		self.check_Argument(sys.argv[1:], "correct number")
		# Non-Mandatory Inputs in Default State:
		output_File_Directory = "./Output_Files/"
		project_Name = "SNV_Project"
		SNV_Title = "'SNV Analysis'"
		low_Storage_Space = False
		profile_SCVs = False
		profile_All = False
		profile_Type = "NT"
		threads = "1"
		bin_Name = "Unknown_Bin_Name" # This Wil Error During Part 2 (Review Code Output to Find Bin Name Options)
		# Looping through arguments given, taking every other argument after Run_Full_Pangenome.py
		print("We Processed the Following Inputs For You: ", sys.argv[1:])
		arg_index = 1
		for _ in sys.argv[1:]:
			arg = sys.argv[arg_index]
			# Double Dash ('--') Means No Argument Follows
			if arg.startswith("--"):
				arg_index += 1
			# Single Dash ('-') Means an Argument Follows (Get The Argument)
			elif arg.startswith("-"):
				user_Input = sys.argv[arg_index + 1]
				arg_index += 2
				self.check_Argument(user_Input, "general")
			# Match Argument to Variable
			if arg == "-in":
				self.check_Argument(user_Input, "file_Found")
				input_Fasta_File_Directory = self.check_Argument(user_Input, "File")
			elif arg == "-out":
				output_File_Directory = self.check_Argument(user_Input, "File")
			elif arg == "-S":
				self.check_Argument(user_Input, "file_Found")
				SNV_Fasta_File_of_Interest = self.check_Argument(user_Input, "Fasta")
			elif arg == "-bin":
				bin_Name = self.check_Argument(user_Input, "bin")
			elif arg == "-name":
				project_Name = self.check_Argument(user_Input, "filename")
			elif arg == "-title":
				SNV_Title = self.check_Argument(user_Input, "Title")
			elif arg == "-T":
				threads = self.check_Argument(user_Input, "int")
			elif arg == "-type":
				profile_Type = self.check_Argument(user_Input, "profile")
			elif arg == "--profile_SCVs":
				profile_SCVs = True
			elif arg == "--low_storage":
				low_Storage_Space = True
			elif arg == "--full_run":
				profile_All = True
			else:
				print("Incorrect Input Prefix. Please Reformat to the Style Above")
				exit()
			if arg_index == len(sys.argv):
				break
		# Tell User If They Didnt Give a Bin
		if bin_Name == "Unknown_Bin_Name":
			print("No Bin Given. You Will Pay Later")
		# Record the User Inputs and Excecute Code if Sufficient Inputs Given
		try:
			return input_Fasta_File_Directory, output_File_Directory, SNV_Fasta_File_of_Interest, SNV_Title, bin_Name, project_Name, threads, profile_Type, profile_SCVs, profile_All, low_Storage_Space
		# Raise Error if One or More of These Fine Variables is Not Defined
		except:
			print("Error: Missing Input File (-in or -S) Flag, which is Required when Running the Program via the Command line")
			print("Mandatory Minimum Expected: $ python Run_Full_Pangenome.py -in [Input File Directory] -S [Fasta File of Interest]")
			exit()
	
	def check_Argument(self, user_Input, stage):
		# Should have an Even Amount of Inputs After the .py File
		if stage == "correct number":
			counting_Input = sum([not need_Pair.startswith("--") for need_Pair in user_Input])
			if counting_Input%2 != 0:
				print("You have an ODD Number of Flag Inputs. Did you Add a Single Dash Flag Without an Input??")
				exit()
		# In General, we do NOT want inputs that start with dashes
		elif stage == "general":
			if user_Input.startswith("-"):
				print("Your Input ", user_Input," Started with a '-', which is NOT Allowed. Please Fix This")
				print("Possible User Error: Did You Place 2 Flags Next to Each Other WITHOUT Their Arguments")
				exit()
		# All Files must end with a Slash
		elif stage == "File":
			if user_Input.endswith("/"):
				return user_Input
			else:
				return user_Input + "/"
		# Check if User Input File Exists
		elif stage == "file_Found":
			if not os.path.exists(user_Input):
				print("Your Input Was: ", user_Input, " And it Does Not Exist")
				exit()
		# Check if User Input File is a Fasta File
		elif stage == "Fasta":
			possible_Fasta_Files = tuple([".fasta",  ".fna", ".fsa", ".mpfa", ".fa"])
			if not user_Input.endswith(possible_Fasta_Files):
				print("Your Input to '-S [Fasta File of Interest]' Was: ", user_Input, " And it is NOT a Fasta File")
				exit()
			return user_Input
		# Check if Bin is There
		elif stage == "bin":
			print("I Have No Way to Check Your Bin Input at This Time. If It Fails, Just Redo the 'SNV_Analysis' Step in Part 2")
			print("If You Are Clueless on What a Bin is, it is the New Contig Name AFTER running 'anvi-script-reformat-fasta' in the fasta_Files Output Folder")
			print("To Check Good Bins Later, Run the 'anvi-import-collection' Step in the Code (Or Just Check The Output of This Run)")
			return user_Input
		# There Should NOT be Spaces in a Filename for the love of coding (Also no bad characters)
		elif stage == "filename":
			for bad_Char in [" ", ";", ")", "(", "-", ".", "?", "/"]:
				user_Input.replace(bad_Char, "_") # Replace with Underscore
			return user_Input.replace(" ", "_") 
		# If There Are Spaces in the Title, It MUST be Surrounded by Quotes
		elif stage == "Title":
			if " " in user_Input:
				if not user_Input.endswith("'"):
					user_Input = user_Input + "'"
				if not user_Input.startswith("'"):
					user_Input = "'" + user_Input
			return user_Input
		# Threads must be an INTEGER
		elif stage == "int":
			if user_Input.isnumeric():
				return user_Input
			else:
				print("You Inputed a Thread That is NOT AN INTEGER. Lol what does that even mean")
				exit()
		# Only Three Options for Profile
		elif stage == "profile":
			if user_Input in ["NT", "AA", "CDN"]:
				return user_Input
			print("The Profile Type Can Only be a Nucleotide: NT, Amino Acid: AA, or Codon: CDN.")
			print("Please Use one of the Following Options: NT, AA, or CDN")
			exit() 
		else:
			print("What Argument am I checking?")
			exit()
		

# ----------------------------------------------------------------------------------------------------------------#

class data_Preparation(helper_Files):
	def __init__(self, input_fasta_file_directory, output_fasta_file_directory):
		super().__init__(output_fasta_file_directory)  # Get Variables Inherited from the helper_Files Class
		self.input_fasta_file_directory = input_fasta_file_directory
		self.dont_touch = set() # If you are nervous, add files you DEFINATELY dont want changed in the variable 'dont_touch'

	def remove_Extended_Nucleotides(self):
		"""
		Removing Non-ATCGN Nucleotides. It should be DNA, so not doing anything with U (unsure if RNA works)
		"""
		print("\nRemoving Non-ATCGN Nucleotides")
		# loop through the files in the directory
		for fasta_file in os.listdir(self.input_fasta_file_directory):
            # Make sure I am only changing FASTA files
			possible_Fasta_Files = tuple([".fasta",  ".fna", ".fsa", ".mpfa", ".fa"])
			if fasta_file not in self.dont_touch and fasta_file.endswith(possible_Fasta_Files):
				# Get the previous name of the File + data
				old_file = self.input_fasta_file_directory + fasta_file
				fasta_sequences = SeqIO.parse(open(old_file), 'fasta')
				# Create new file and write to it
				base = os.path.splitext(fasta_file)[0]
				new_file = self.fasta_Files + base + ".fasta"
				
				records = []
				with open(new_file, "w+") as out_file:
					# Modify the contigs one by one (Using enumerate as online it said it needed to be a generator ??)
					for record in fasta_sequences:
						#print("old \n", record)
						# Get Sequence of the Contig
						new_sequence = str(record.seq)
						# Removing Extented Nucleotides Found via https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2865858/
						for NT in ["R", "Y", "M", "K", "S", "W", "H", "B", "V", "D"]:
							new_sequence = new_sequence.replace(NT,"N")
						record.seq = Seq(new_sequence, record.seq.alphabet)
						#print("new \n", old_record) 
						records.append(record)
						# Write the File
					SeqIO.write(records, new_file, "fasta")
			else:
				print("Not Using: ", fasta_file)

	# Fix File Names Loop
	def fix_file_names(self):
		"""
		ANVIO Databases only uses files with a '.fa' extension, NO spaces/non-ASCII characters, and where the first digit is NOT a number
		The following function will check for this and make appropriate changes to the FASTA file names
		"""
		print("\nFixing Names of Files")
		# loop through the files in the directory
		for fasta_file in list(os.listdir(self.fasta_Files)):
			# Make sure I am only changing FASTA files
			fixable_ext = tuple([".fasta",  ".fna", ".fsa", ".mpfa", ".fa"])
			if fasta_file not in self.dont_touch and fasta_file.endswith(fixable_ext):
				# Edit the File Names to Work with ANVIO Code
				old_file = self.fasta_Files + fasta_file
				new_file = self.filename_Edits(fasta_file)
				print(old_file + "   ->   " + new_file)
				os.rename(old_file, new_file)
			else:
				print("Not Using: ", fasta_file)
	
	# Fix One File's Names
	def filename_Edits(self, fasta_file):
		# Get the previous name of the File
		old_file = self.fasta_Files + fasta_file
		base = os.path.splitext(fasta_file)[0]
		# Please limit the characters to ASCII letters, digits (not in first position), and the underscore ('_') character
		# This ALSO means no spaces
		for bad_Char in ["'"]:
			base = base.replace(bad_Char, "")  # Delete Character
		for bad_Char in [" ", ";", ")", "(", "-", ".", "?", "/"]:
			base = base.replace(bad_Char, "_") # Replace with Underscore
		# If First letter is a digit, add string
		if base[0].isdigit():
			base = "ANVIO_" + base
		# Change File Extension, must be .fa for ANVIO
		new_file = self.fasta_Files + base + '.fa'
		return new_file


# ----------------------------------------------------------------------------------------------------------------#

class ANVIO_Storage(helper_Files):
	# Initialize Files
	def __init__(self, output_File_Directory, project_name, threads):
		super().__init__(output_File_Directory)  # Get Variables Inherited from the helper_Files Class
		self.dont_touch = {}  # If you are nervous, add files you DEFINATELY dont want changed in the variable 'dont_touch'
		self.project_name = project_name
		self.threads = int(threads)

	def storage_Step(self, stage):
		# Divide up the work with threading. Threads per loop will be optimally decided to maximize parrell processes
		max_Useful_Pool = sum([Fasta_File.endswith(".fa") for Fasta_File in os.listdir(self.fasta_Files)]) # If the user has more threads than files, its not worth wasting them
		if stage == "format":
			max_Useful_Pool = 1 # Will Decide Later if It Truely is 'Better' as it is a fast step anyways
		threads_per_loop = min(16, max(self.threads//max_Useful_Pool, 1))  # Not worth going over 16 as the extra threads are not useful
		best_Pool = min(max_Useful_Pool, self.threads)

		pool = multiprocessing.Pool(best_Pool)  # The number of parrelel loops. It helps to maximize this. Note: this will make a lot of intermediate storage
		# Looping through Directory to get Fasta Files (making list() to ensure we use only files PREVIOUSLY present)
		for filename in list(os.listdir(self.fasta_Files)):
			# Make sure I am only using the correct files (NOTE: '-clean.fa' will appear and disappear as the program goes)
			if filename not in self.dont_touch and filename.endswith(".fa") and not filename.endswith("-clean.fa"):
				# Threading means that all loops will be executed at once (cannot rely on each other)
				if stage == "format":
					pool.apply_async(self.format_Fasta_Files, args=(filename,))
				elif stage == "convert":
					pool.apply_async(self.create_Contig_Database, args=(filename, str(threads_per_loop)))
				else:
					print("This Stage Does Not Exist")
					exit()
			else:
				print("Not ", stage, "ing: ", filename)
		pool.close()
		pool.join()

	def format_Fasta_Files(self, fasta_Filename):
		"""
		The following function will remove contigs that have big gaps or small length. It can also simplify the contig names.
		"""
		print("\nReformating Fasta")
		# NOTE: After 2 hours of crying, I finally realized that you will get a 'Not a Fasta' error if you try to overwrite the file (input_File != outfile)
			# The error: anvio.fastalib.FastaLibError: Fasta Lib Error: File './Output_Files/JC_0293_contigs.fa' does not seem to be a FASTA file.
		# Add the '--min-len' flag if you want to exlude contigs below a certain length. Default is 0
		# Add the '--max-percentage-gaps' flag if you want to exlude contigs that have a gap above a certain percent. Default is 100
		# Add the '--simplify-names' flag if you want to simplify contig/scaffold names ('scaffold 1' -> '001'). It helps as the BAM converter will do this with bad deflines arbitrarily
		# Add the '--report-file fasta_deflines_map.txt' flag if you used '--simplify-names' and want to map what was changed (basically use it if you add '--simplify-names')
		min_length = "0"  # Default is 0
		max_gaps = "100"     # Default is 100 (%)
		base = os.path.splitext(fasta_Filename)[0]
		defline_map = self.contigs_Database + base + ".txt"
		input_File = self.fasta_Files + fasta_Filename
		output_File = self.fasta_Files + "hold_" + fasta_Filename
		command = "anvi-script-reformat-fasta '" + input_File + "' --min-len " + min_length + " --max-percentage-gaps " + max_gaps + " -o " + output_File + " --simplify-names --report-file " + defline_map
		self.execute_Command(command)
		# Only keep one of the files
		command = "mv " + output_File + " " + input_File
		self.execute_Command(command)

	def create_Contig_Database(self, filename, threads_per_loop):
		print("\nConverting FASTA to ANVIO Database")
        # Convert FASTA format to ANVIO Database format
		# Add the '--skip-mindful-splitting' flag to not search for optimal splitting and just do the default (makes it run faster)
		# Add the '--split-length' flag to specify which legnth to split (ANVIO will still look for best if '--skip-mindful-splitting' flag is not present. Default is 20,000
		base = os.path.splitext(filename)[0]
		database = base + ".db"
		filename = self.fasta_Files + filename
		command = "anvi-gen-contigs-database -f " + filename + " -o " + self.contigs_Database + database + " -n " + self.project_name
		self.execute_Command(command)
		# Run Hmms
		command = "anvi-run-hmms -c " + self.contigs_Database + database + " --num-threads " + threads_per_loop
		self.execute_Command(command)
        # Annotate the genome: optional
		command = "anvi-run-ncbi-cogs -c " + self.contigs_Database + database + " --num-threads " + threads_per_loop  # This requires you to have previously installed 'anvi-setup-ncbi-cogs'
		self.execute_Command(command)


# ----------------------------------------------------------------------------------------------------------------#

class Quality_Filtering(helper_Files):
	# Based on the ANVIO outline: http://merenlab.org/tutorials/assembly-based-metagenomics/#quality-filtering
	def __init__(self, input_Reads_Directory, output_File_Directory, threads):
		super().__init__(output_File_Directory)  # Get Variables Inherited from the helper_Files Class
		self.input_Reads_Directory = input_Reads_Directory
		self.threads = threads
		self.dont_use = {}

	def unzip_Reads(self):
		# Cannot Filter Zipped Files (I think?)
		print("\nUnzipping Read Files")
		command = "gunzip " + self.input_Reads_Directory + "*.fastq.gz"
		self.execute_Command(command)

	def store_Reads(self):
		# Make a tab delimited file containing the read file paths that we are filtering
		print("\nCreating Reads Database")
		with open(self.filter_Directory + "samples.txt", 'w+') as Reads_File:
            # Make file tab delimited, include ANVIO's header
			writer = csv.writer(Reads_File, delimiter='\t')
			writer.writerow(["sample", "r1", "r2"])

			# loop through the files in the directory to change them
			for Read in os.listdir(self.input_Reads_Directory):
				# Make sure I am only changing the correct files
				if Read not in self.dont_use and Read.endswith(".fastq") and "R1" in Read:
					# Find Sample Name and Forward/Reverse Reads
					Split_Read = Read.split("R1")
					R1 = self.input_Reads_Directory + Split_Read[0] + "R1" + Split_Read[1]
					R2 = self.input_Reads_Directory + Split_Read[0] + "R2" + Split_Read[1]
					sample_name = Split_Read[0]
					# Write the File
					writer.writerow([sample_name, R1, R2])
				else:
					print("Not Using: ", Read)
	
	def filter_Reads(self):
		# Create config files for illumina-utils
		command = "iu-gen-configs " + self.filter_Directory + "samples.txt -o " + self.filter_Directory
		self.execute_Command(command)
		
		# Divide up the work with threading. Threads per loop will be optimally decided to maximize parrell processes
		max_Useful_Pool = sum([filter_File.endswith(".ini") for filter_File in os.listdir(self.filter_Directory)]) # If the user has more threads than files, its not worth wasting them
		best_Pool = min(max_Useful_Pool, int(self.threads))

		pool = multiprocessing.Pool(best_Pool)  # The number of parrelel loops. It helps to maximize this. Note: this will make a lot of intermediate storage
		# Looping through Directory to get Fasta Files (making list() to ensure we use only files PREVIOUSLY present)
		for filename in os.listdir(self.filter_Directory):
			# Make sure I am only filtering the correct files
			if filename not in self.dont_use and filename.endswith(".ini"):
				# Threading means that all loops will be executed at once (cannot rely on each other)
				pool.apply_async(self.apply_Filter, args=(filename,))
			else:
				print("Not Filtering: ", filename)
		pool.close()
		pool.join()

	def apply_Filter(self, filename):
		# Filter the Reads
		# If using overlapping reads, use 'iu-merge-pairs' instead
		command = "iu-filter-quality-minoche " + self.filter_Directory + filename
		self.execute_Command(command)


# ----------------------------------------------------------------------------------------------------------------#

class add_BAM(helper_Files):
	def __init__(self, output_File_Directory, Fasta_File, threads, profile_SCVs, profile_All, low_Storage_Space = True):
		super().__init__(output_File_Directory)  # Get Variables Inherited from the helper_Files Class
		self.profile_SCVs = profile_SCVs
		self.threads = int(threads)
		self.low_Storage_Space = low_Storage_Space
		self.profile_One_Fasta_File = Fasta_File
		self.profile_All = profile_All
		input_Files_For_SNV, output_File_Directory, Fasta_File_For_SNV, threads, profile_SCVs, profile_All, low_Storage_Space

	def task(self, stage):
		print("\nGenerating the BAM Files")
		dont_use = {}
		# Get samples
		with open(self.filter_Directory + "samples.txt") as sample:
			reader = csv.reader(sample, delimiter="\t")
			samples = list(zip(*reader))[0]
		samples = samples[1:] # remove the header 'sample'

		# Divide up the work with threading. Threads per loop will be optimally decided to maximize parrell processes
		max_Useful_Pool = sum([Fasta_File.endswith(".fa") for Fasta_File in os.listdir(self.fasta_Files)]) # If the user has more threads than files, its not worth wasting them
		if self.profile_All == False:
			max_Useful_Pool = 1
		threads_per_loop = min(32, max(self.threads//max_Useful_Pool, 1)) # Max useful threads for this process is 16 (really 8)
		best_Pool = min(max_Useful_Pool, self.threads)
		if self.low_Storage_Space and stage == "generate":
			best_Pool = 1  # As the 'generate' process with Bamtools will create ~12G of extra storage per loop (to be deleted at end), you may not have the extra room for parellel processing
			threads_per_loop = min(self.threads, 16)  # Max useful threads for this process is 16 (really 8)

		pool = multiprocessing.Pool(best_Pool)  # The number of parrelel loops. It helps to maximize this. Note: this will make a lot of intermediate storage
		# loop through the files in the directory to change them
		for Fasta_File in os.listdir(self.fasta_Files):
			# If User Only Wants to profile One Fasta File
			print(self.profile_One_Fasta_File, " try: ", Fasta_File)
			if self.profile_All == False and self.profile_One_Fasta_File != Fasta_File:
				continue
			# Make sure I am only changing the correct files
			if Fasta_File not in dont_use and Fasta_File.endswith(".fa"):
				# Threading means that all loops will be executed at once (cannot rely on each other)
				if stage == "generate":
					pool.apply_async(self.generate, args=(Fasta_File, str(threads_per_loop), samples))
				elif stage == "anvi_sort_index":
					pool.apply_async(self.anvi_Sort_Index_BAM, args=(Fasta_File, str(threads_per_loop), samples))
				elif stage == "profile":
					pool.apply_async(self.profile, args=(Fasta_File, str(threads_per_loop), samples))
				elif stage == "merge":
					pool.apply_async(self.merge, args=(Fasta_File,))
			else:
				print("Not Using: ", Fasta_File)   
		pool.close()
		pool.join()

	def generate(self, Fasta_File, generate_threads, samples):
		"""
		Uses Bowtie2 to map '.fa' to '.sam', SAMTools to map '.sam' to '.bam'
		"""
		# Make an Output Folder for This Bowtie
		bowtie_Index = os.path.splitext(Fasta_File)[0]
		bowtie_Folder = self.bam_Directory + bowtie_Index + "/"
		command = "mkdir " + bowtie_Folder
		self.execute_Command(command)
		# Buid a Bowtie2 Database
		# Add the '--quiet' flag to stop a hundred lines from printing out (remove if debugging)
		command = "bowtie2-build --quiet " + self.fasta_Files + Fasta_File + " " + bowtie_Folder + bowtie_Index
		self.execute_Command(command)
		

		# loop through the files in the directory to change them
		dont_use = {'JC-2272_S289_L004_','JC-2285_S258_L004_','JC-2289_S301_L004_'}
		for sample_name in samples:
			if sample_name in dont_use:
				continue
			R1 = sample_name + "-QUALITY_PASSED_R1.fastq"
			R2 = sample_name + "-QUALITY_PASSED_R2.fastq"
			base = os.path.splitext(sample_name)[0]
			print("Looping through samples:\n", sample_name, "\n", R1, "\n", R2)

			# do the bowtie mapping to get the SAM file:
			# it uses '--phred33' by default but you can do '--phred64' (I believe most data now recommends 33 over 64)
			# Add the '-q' flag to specify that input files are '.fastq'
			command = "bowtie2 -q --phred33 --threads " + generate_threads + " -x " + bowtie_Folder + bowtie_Index + " -1 " + self.filter_Directory + R1 + " -2 " + self.filter_Directory + R2 + " --no-unal -S " + bowtie_Folder + base + ".sam"
			self.execute_Command(command)

			# covert the resulting SAM file to a BAM file:
			command = "samtools view --threads " + generate_threads + " -F 4 -bS " + bowtie_Folder + base + ".sam > " + bowtie_Folder + base + "-RAW.bam"
			self.execute_Command(command)

			# Sort BAM File:
			# NOTE: Something like '[bam_sort_core] merging from 8 files and 1 in-memory blocks...' is not an error (at least I have found it not to be)
			command = "samtools sort --threads " + generate_threads + " " + bowtie_Folder + base + "-RAW.bam -o " + bowtie_Folder + base + ".bam"
			self.execute_Command(command)
			
			# Index BAM File
			# Add the flag '-@' to specify thread count (idk why its an @ symbol lol)
			command = "samtools index " + bowtie_Folder + base + ".bam -@ " + generate_threads
			self.execute_Command(command)

			# Remove Intermediate Files (as they are BIG):
			command = "rm -r " + bowtie_Folder + base + ".sam " + bowtie_Folder + base + "-RAW.bam"
			self.execute_Command(command)
			
		# Delete Extra Bowtie2 Files Created in Generation Stage
		command = "rm -r " + self.bowtie_Folder + "*.bt2"
		self.execute_Command(command)
			
	def anvi_Sort_Index_BAM(self, Fasta_File, profile_threads, samples):
		print("\nProfiling BAMs")
		# Profile the database
		bowtie_Index = os.path.splitext(Fasta_File)[0]
		bowtie_Folder = self.bam_Directory + bowtie_Index + "/"

		# loop through the files in the directory to change them
		for sample_name in samples:
			base = os.path.splitext(sample_name)[0]
			print("Looping through samples:\n", sample_name)

			# Process BAM files for ANVIO's Database
			# Add the '--profile-SCVs' flag if you want them profiled
			# Add the '--skip-SNV-profiling' if you dont care for SNVs
			# Add the '-M #' flag where the # is the minimum length of a contig to analyze in your BAM (default is 1000)
			command = "anvi-init-bam -T " + profile_threads + " " + bowtie_Folder + base + ".bam -o " + bowtie_Folder + base + "_new.bam"
			self.execute_Command(command)

	def profile(self, Fasta_File, profile_threads, samples):
		print("\nProfiling BAMs")
		# Profile the database
		genome_Isolate_Base = os.path.splitext(Fasta_File)[0]  # Also the Bowtie Index (Meaning ~Different But Same)
		bowtie_Index = os.path.splitext(Fasta_File)[0]
		bowtie_Folder = self.bam_Directory + bowtie_Index + "/"
		add_SCVs = ""
		if self.profile_SCVs:
			add_SCVs = " --profile-SCVs"

		# loop through the files in the directory to change them
		for sample_name in samples:
			base = os.path.splitext(sample_name)[0]
			print("Looping through samples:\n", sample_name)

			# Process BAM files for ANVIO's Database
			# Add the '--profile-SCVs' flag if you want them profiled
			# Add the '--skip-SNV-profiling' if you dont care for SNVs
			# Add the '-M #' flag where the # is the minimum length of a contig to analyze in your BAM (default is 1000)
			min_Contig_length = "1000"
			command = "anvi-profile -c " + self.contigs_Database + genome_Isolate_Base + ".db -i " + bowtie_Folder + base + ".bam -M " + min_Contig_length + add_SCVs + " --num-threads " + profile_threads + " -o " + bowtie_Folder + base
			self.execute_Command(command)

	def merge(self, Fasta_File):
		print("\Merging BAM Profiles")
		# Merge the databases
		genome_Isolate_Base = os.path.splitext(Fasta_File)[0]  # Also the Bowtie Index (Meaning ~Different But Same)
		bowtie_Index = os.path.splitext(Fasta_File)[0]
		bowtie_Folder = self.bam_Directory + bowtie_Index + "/"

		# Done on ALL the profiles
		command = "anvi-merge " + bowtie_Folder + "*/PROFILE.db -o " + bowtie_Folder + genome_Isolate_Base + "-MERGED -c " + self.contigs_Database + genome_Isolate_Base + ".db"
		self.execute_Command(command)


# ----------------------------------------------------------------------------------------------------------------#

class supply_Collection(helper_Files):
	def __init__(self, output_File_Directory, Fasta_File, profile_All):
		super().__init__(output_File_Directory)  # Get Variables Inherited from the helper_Files Class
		self.profile_One_Fasta_File = Fasta_File
		self.profile_All = profile_All

	def generate_Collection(self):
		dont_use = {}
		# loop through the files in the directory to change them
		for Fasta_File in os.listdir(self.fasta_Files):
			# If User Only Wants to profile One Fasta File
			if self.profile_All == False and self.profile_One_Fasta_File != Fasta_File:
				continue
			# Make sure I am only changing the correct files
			if Fasta_File not in dont_use and Fasta_File.endswith(".fa"):
				# Get Contig Database that was Created from Isolate During Storage
				genome_Isolate_Base = os.path.splitext(Fasta_File)[0]  # coinciently also the Bowtie Index (for ease)
				contig_Database = self.contigs_Database + genome_Isolate_Base + ".db"
				
				with open(self.collection_Files + genome_Isolate_Base + "-GENOME-COLLECTION.txt", 'w+') as collection_File:
				    # Make file tab delimited
					writer = csv.writer(collection_File, delimiter='\t')

					# Use SQLite3 to loop over contig
					conn = sqlite3.connect(contig_Database)
					for split_name in conn.execute('select split from splits_basic_info'):
						split_name = split_name[0]
						genome_Name = split_name.split("_split_")[0]
						# print(split_name, "\t", genome_Name)
						writer.writerow([split_name, genome_Name])
						conn.commit()

	def import_Collection(self):
		dont_use = {}
		# loop through the files in the directory to change them
		for Fasta_File in os.listdir(self.fasta_Files):
			# If User Only Wants to profile One Fasta File
			if self.profile_All == False and self.profile_One_Fasta_File != Fasta_File:
				continue
			# Make sure I am only changing the correct files
			if Fasta_File not in dont_use and Fasta_File.endswith(".fa"):
				# Get Contig Database that was Created from Isolate During Storage
				genome_Isolate_Base = os.path.splitext(Fasta_File)[0]  # coinciently also the Bowtie Index (for ease)
				contig_Database = self.contigs_Database + genome_Isolate_Base + ".db"
				collection_File = self.collection_Files + genome_Isolate_Base + "-GENOME-COLLECTION.txt"
				collection_Name = "Genomes"
				bowtie_Index = os.path.splitext(Fasta_File)[0]
				bowtie_Folder = self.bam_Directory + bowtie_Index + "/"
				merged_Profile_Database = bowtie_Folder + genome_Isolate_Base + "-MERGED/PROFILE.db"
				collection_Output = self.collection_Files + "collections_and_bins.txt"
				
				# Import the Collection
				command = "anvi-import-collection " + collection_File + " -c " + contig_Database + " -p " + merged_Profile_Database + " -C " + collection_Name
				self.execute_Command(command)

				# Check to see if it worked
				command = "nohup anvi-show-collections-and-bins -p " + merged_Profile_Database + " > " + collection_Output
				self.execute_Command(command)

				# Summarize Results. Not Really Needed
				summary_Output = self.collection_Files + genome_Isolate_Base + "-SUMMARY"
				command = "anvi-summarize -c " + contig_Database + " -p " + merged_Profile_Database + " -C " + collection_Name + " --init-gene-coverages -o " + summary_Output
				self.execute_Command(command)


# ----------------------------------------------------------------------------------------------------------------#

class create_SNV(helper_Files):
	def __init__(self, output_File_Directory, Fasta_File, profile_All):
		super().__init__(output_File_Directory)  # Get Variables Inherited from the helper_Files Class
		self.profile_One_Fasta_File = Fasta_File
		self.profile_All = profile_All

	def variability_Profile(self, profile_Type, bin_Name):
		print("\nCreating Variability Profile")
		dont_use = {}
		# loop through the files in the directory to change them
		for Fasta_File in os.listdir(self.fasta_Files):
			print(self.profile_One_Fasta_File, "     ", Fasta_File)
			if self.profile_All == False and self.profile_One_Fasta_File != Fasta_File:
				continue
			# Make sure I am only changing the correct files
			if Fasta_File not in dont_use and Fasta_File.endswith(".fa"):
				# Get Contig Database that was Created from Isolate During Storage
				genome_Isolate_Base = os.path.splitext(Fasta_File)[0]  # coinciently also the Bowtie Index (for ease)
				contig_Database = self.contigs_Database + genome_Isolate_Base + ".db"
				collection_Name = "Genomes"
				bowtie_Index = os.path.splitext(Fasta_File)[0]
				bowtie_Folder = self.bam_Directory + bowtie_Index + "/"
				merged_Profile_Database = bowtie_Folder + genome_Isolate_Base + "-MERGED/PROFILE.db"

				min_Coverage = "20"  # Minimum Coverage in Genomes
				min_Occurrence = "1" # Minimum occurence in genomes (3 seems like a good number)
				command = "anvi-gen-variability-profile -c " + contig_Database + " -p " + merged_Profile_Database + " --engine " + profile_Type + " --min-coverage-in-each-sample " + min_Coverage + " --min-occurrence " + min_Occurrence + " --include-split-names --quince-mode -C " + collection_Name + " -b " + bin_Name + " -o " + self.output_SNV + genome_Isolate_Base + "_SNV.txt"
				self.execute_Command(command)

	def gen_SNV(self):
		print("\nGenerating SNV")
		dont_use = {}
		# loop through the files in the directory to change them
		for Fasta_File in os.listdir(self.fasta_Files):
			if self.profile_All == False and self.profile_One_Fasta_File != Fasta_File:
				continue
			# Make sure I am only changing the correct files
			if Fasta_File not in dont_use and Fasta_File.endswith(".fa"):
				# Get Contig Database that was Created from Isolate During Storage
				genome_Isolate_Base = os.path.splitext(Fasta_File)[0]  # coinciently also the Bowtie Index (for ease)

				command = "anvi-script-snvs-to-interactive " + self.output_SNV + genome_Isolate_Base + "_SNV.txt -o " + self.output_SNV + genome_Isolate_Base + "/"
				self.execute_Command(command)

	def display_SNV(self, SNV_Title, Fasta_File):
		genome_Isolate_Base = os.path.splitext(Fasta_File)[0]  # coinciently also the Bowtie Index (for ease)
		output_SNV_Folder = self.output_SNV + genome_Isolate_Base + "/"
		command = "sudo anvi-interactive --manual -p " + output_SNV_Folder + "profile.db --tree " + output_SNV_Folder + "tree.txt --view-data " + output_SNV_Folder + "view.txt --title " + SNV_Title + " --manual"
		self.execute_Command(command)


# ----------------------------------------------------------------------------------------------------------------#

if __name__ == "__main__":
	"""
	Notes of Using this File:
		1. input_Files_For_SNV: This folder should have all the fasta and reads files of the MAGs inside. The Reads and Fasta Files must all start with JC-# (really *-*).
								If this is not the case, you will need to change the create_BAM class to find a way to match up fasta files to R1/R2 reads.
								The Fasta Files can end with ".fasta",  ".fna", ".fsa", ".mpfa", ".fa" and will later be changed to '.fa'.
								The read files can end with either '.fastq.gz' (where it will be unziped) or .'fastq'.
								If it is not obvious, Read files CANNOT have spaces in them (I handle the case if FASTA files have spaces)
		2. output_File_Directory: All '.fa' files come from the data_Preparation class
								  All of the contigs_Database come from the ANVIO_Storage class
								  All of the Quality_Filtering_Files files come from the Quality_Filtering class
								  All of the BAM_Files come from the create_BAM class
		3. with_BAM.task("generate"): If you are stuck on one of the Samtools step (for >2 hours), it is probably because you dont have enough memory (program will just pause without error)
	"""
	# Gather User Inputs (From Command Line)
	get_Info = user_Inputs()
	input_Files_For_SNV, output_File_Directory, Fasta_File_For_SNV, SNV_Title, bin_Name, project_Name, threads, profile_Type, profile_SCVs, profile_All, low_Storage_Space = get_Info.gather_Inputs()
	# If No Inputs to Command Line: Use Inputs Below
	if len(sys.argv) == 1:
		# Define Files Your Files Here
		Fasta_File_For_SNV = "Input_Files_BP/SAF32_gene_NT.fa"     # Name of the Fasta File you Want to have an SNV Profile For (To be Compared Against Other Input Files). Leave as "" if Unknown
		output_File_Directory = "./Output_Files_BP/"   		# Must end with '/'. Files could be edited (or not).
		input_Files_For_SNV = "./Input_Files_BP/"      		# Must end with '/'. Contains Fasta Files + Reads 1/2. NOTHING here will be edited
		project_Name = "Bacillus SNV"      			    # Name of Project
		SNV_Title = "'Bacillus SNV'"		           		# Title of the final SNV Profile
		low_Storage_Space = False  				           # Generation of BAM files is a VERY storage-costly process (space that will be deleted at end). Say True if you are limited on Space
		profile_All = False               					# If True, Will Perform 'Part One' on ALL Fasta Inputs, Ignoring Fasta File Chosen (Unneccesary if You Know Which File You Want)
		threads = 5                                		# specify the maximum number of threads you are willing to use on this VERY TIME-CONSUMING program
		# Working on Adding! Do Not Edit!
		profile_Type = "NT"  #KEEP THIS NT (ERROR!)		# The engine can focus on nucleotides (NT), codons (CDN), or an amino acids (AA). If CDN/AA, then profile_SCVs MUST be True
		profile_SCVs = False #KEEP THIS FALSE (ERROR!!)  	# Profiles SCV (needed if you do CDN/AA). Can have just in case you want to profile one day without redoing work in 'add_BAM'
	
	# Make Output Folders in output_File_Directory
	header = helper_Files(output_File_Directory)
	header.make_Output_Folders(output_File_Directory)
	
	
	# -----------------------------------------------PART ONE: No Human Needed-----------------------------------------------------------------#
	Fasta_File_For_SNV =  Fasta_File_For_SNV.split("/")[-1]  # Get the Filename from File's Path
	# Here to capture your mistake if you want AA/CDN but didnt have SCVs as true
	if profile_Type == "AA" or profile_Type == "CDN":
		profile_SCVs = True

	if False:
		# Data Preparation
		edit_FASTA_Files = data_Preparation(input_Files_For_SNV, output_File_Directory)  # Create Class Instance
		edit_FASTA_Files.remove_Extended_Nucleotides() 	# Remove non-ATCGN Nucleotides (ANVIO Cannot Deal With Them)
		edit_FASTA_Files.fix_file_names()    				# Edit Fasta File names to comply with ANVIO database
		Fasta_File_For_SNV = edit_FASTA_Files.filename_Edits(Fasta_File_For_SNV) # Edit Input File for SNV's Name to Match How we Changed the Actual File's Name
		Fasta_File_For_SNV =  Fasta_File_For_SNV.split("/")[-1]  # Get the Filename from File's Path

		# Create Database Storage with ANVIO
		Make_Storage = ANVIO_Storage(output_File_Directory, project_Name, threads)  # Create Class Instance.
		Make_Storage.storage_Step("format")   # Remove Bad Contigs and Format Fasta File to Work with Anvio's Software
		Make_Storage.storage_Step("convert")  # Convert FASTA '.fa' to '.db', YOU MUST MUST run 'anvi-setup-ncbi-cogs' once in the terminal
		
		# Filter the Reads given from Asembly (R1/R2). I have set up the code for paired end reads in mind (check if you are different)
		Quality_Control = Quality_Filtering(input_Files_For_SNV, output_File_Directory, threads)  # Create Class Instance
		Quality_Control.unzip_Reads()   # Unzip the '.fastq.gz' files. If none present, nothing happens
		Quality_Control.store_Reads()   # Stores Files likes '*R1*.fastq' and '*R2*.fastq' for filtering in samples.txt
		Quality_Control.filter_Reads()  # Filters all files in the samples.txt

		# Create Bam File
		with_BAM = add_BAM(output_File_Directory, Fasta_File_For_SNV, threads, profile_SCVs, profile_All, low_Storage_Space)  # Create Class Instance
		with_BAM.task("generate")  		# Generates an indexed and sorted BAM file for each metagenome using the Filtered Reads
		with_BAM.task("profile")  		# Process the BAM files to generate ANVIO profile databases containing the coverage and detection statistics of each genome isolate in a given metagenome
		with_BAM.task("merge")    		# Create merged ANVIO profile database that describes the coverage and detection statistics for each genome isolate across all metagenomes
		
		# Add Collection
		link_Contigs = supply_Collection(output_File_Directory, Fasta_File_For_SNV, profile_All)  # Create Class Instance
		link_Contigs.generate_Collection()  
		link_Contigs.import_Collection()

	# -----------------------------------------------PART TWO: Human Needed-----------------------------------------------------------------#
	bin_Name = "c_000000000033"  # Specify the Bin You Want to See Your SNV Profile On (Can Change Anytime)
	if True:
		# Data Preparation
		edit_FASTA_Files = data_Preparation(input_Files_For_SNV, output_File_Directory)  # Create Class Instance
		Fasta_File_For_SNV = edit_FASTA_Files.filename_Edits(Fasta_File_For_SNV) # Edit Input File for SNV's Name to Match How we Changed the Actual File's Name
		Fasta_File_For_SNV =  Fasta_File_For_SNV.split("/")[-1]  # Get the Filename from File's Path
        
		# Generate SNV
		SNV_Analysis = create_SNV(output_File_Directory, Fasta_File_For_SNV, profile_All)
		#SNV_Analysis.variability_Profile(profile_Type, bin_Name)
		#SNV_Analysis.gen_SNV()
		SNV_Analysis.display_SNV(SNV_Title, Fasta_File_For_SNV)





