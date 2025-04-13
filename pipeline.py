import os 
import subprocess
import random 

# Joshua's Code - Downloading Genomes 
def download_genomes():
    if "Genomes" not in os.getcwd(): #if you aren't already in the Genomes folder
        if not os.path.isdir("Genomes"): #make a directory to store genomes if it doesn't already exist 
            os.system("mkdir Genomes")
        os.chdir("Genomes") #move to the Genomes folder

    accession = ["GCF_014961145.1","GCF_028532485.1","GCF_021391435.1","GCF_004379335.1"] # genome accession codes 

    for x in accession: #loop through accessions
        if x not in ",".join(os.listdir()): #no need to redownload 
            os.system(f"datasets download genome accession {x} --include genome") #download genome from refseq
            os.system("unzip ncbi_dataset.zip") #unzip
            os.system("rm ncbi_dataset.zip") #delete zip file which has generic name to prevent overwrite issues
            os.system("mv ncbi_dataset/data/GCF*/GCF* .") #extract the genome file and move it to current folder, which should be Genomes
            os.system("rm -rf ncbi_dataset/data") #delete the unzipped folder which has unnecessary extra stuff
            os.system("rm README.md") #delete extra file to avoid asking for overwrite
            os.system("rm md5sum.txt") #delete extra file to avoid asking for overwrite
            os.system("rmdir ncbi_dataset") #remove empty folder
    os.chdir("..") #return to main folder

# Hillary's Code - Generating and Inserting Repeats 
# This function is still a work in progress and needs more editting 
def generate_and_insert_repeats():
    if not os.path.isfile("Motifs/ip.txt"):
        bases = ["A", "C", "G", "T"] #nucleotide bases

        # 100 BP random repetitive element 
        motif1 = random.choices(bases, k=100) #100bp random sequence
        seq = "".join(motif1) #join the bases together to make a string

        # 500 BP random repetitive element
        motif2 = random.choices(bases, k=500) #500bp random sequence
        seq2 = "".join(motif2) #join the bases together to make a string

        if not os.path.isdir("Motifs"): #make a directory to store the motifs if it doesn't already exist 
            os.system("mkdir Motifs")
            os.chdir("Motifs") #move to that directory
            # Write the repetitive elements to files - might be necessary/helpful to see the sequences for future analysis 
            with open("motif1.txt", "w") as f:
                f.write(seq)
            with open("motif2.txt", "w") as f:
                f.write(seq2)
            os.chdir("..") #move back to the original directory

        
        # Create a directory to store the modified genomes 
        modified_dir = "Modified_Genomes"
        if not os.path.isdir(modified_dir): #make a directory to store modified genomes if it doesn't already exist
            os.system(f"mkdir {modified_dir} ") #create directory if it doesn't exist
        
        # Create a dictionary of motifs
        motifs = {"motif1": seq, "motif2": seq2} 
        
        num_insertions = [2, 3, 4, 5] #number of insertions to make

        ip = [] #list of insertion points
        mingenomelength = 1000000000000 #1 trillion, arbitrary maximum
        for m in os.listdir("Genomes"): #loop through files
            with open("Genomes/{}".format(m), "r") as f: #open files
                dat = f.read() #get length
                genomelength = len(dat)
                if genomelength < mingenomelength: #get shortest genome
                    mingenomelength = genomelength
        for n in range(len(num_insertions)+1): #insert based on number of insertions
            ip.append(random.randint(0,mingenomelength)) #generate numbers that don't exceed the shortest genome
        
        ip = sorted(ip) #sort the list of insertion points
        print(ip) #print the list of insertion points

        # Write the insertion points to a file - this is optional but can be helpful for debugging (able to check exact positions where the motifs were inserted)
        with open("Motifs/ip.txt", "w") as f:
            for i in ip:
                f.write(f"{i}\n")

        for file in os.listdir("Genomes"): #loop through the genomes
            with open(f"Genomes/{file}", "r") as f:
                genome = f.read()

            accession = file.split(".")[0]

            for motif_name, sequence in motifs.items():
                for count in num_insertions: #loop through the number of insertions
                    mod_genome = genome
                    for i in range(count):
                        insertion_point = ip[i] #iterate through the predetermined list of random insertion points based on the number of repeats
                        mod_genome = mod_genome[:insertion_point] + sequence + mod_genome[insertion_point:]
                    output_filename = f"{accession}_{motif_name}_{count}.fna"
                    
                    with open(f"{modified_dir}/{output_filename}", "w") as out_f:
                        out_f.write(mod_genome)



# Joshua's Code - Running ART (Artificially Simulated Genomes)
def run_art():
    #modified genomes list
    genomes = []
    for file in os.listdir("Modified_Genomes"):
        genomes.append(file)

    if not os.path.isdir("artgens"): #directory for simulated genomes
        os.system("mkdir artgens") #create directory if it doesn't exist
    os.chdir("artgens") #move to that directory

    for gen in genomes:
        #for depth in range(10,110,10): #increase depth from 10 to 100 in increments of 10
        for depth in [10,100]: #runs art with depth of 10 and 100
            out = str(gen.split(".")[0]) + "_" + str(depth) + "_"
            #art flags: -p = paired ends, -na = don't output alignment file, -i = input, -l = read length, -f = fold coverage, -o = output prefix, -m = mean fragment length, -s = standard deviation
            if not os.path.isfile("{}1.fq".format(out)):
                os.system("art_illumina -p -na -i ../Modified_Genomes/{0} -l 151 -m 200 -s 10 -f {1} -o {2}".format(gen,depth,out))
    os.chdir("..") #move back to the original directory
    

# Kimia's Code - Running SPAdes and Unicycler
def run_spades():
    if not os.path.isdir("Spades_Output"):
        #directory to input files
        input_dir = "artgens"

        # SPAdes Run 
        #making output file
        if not os.path.isdir("Spades_Output"): #make a directory to store SPAdes output if it doesn't already exist 
            os.system("mkdir Spades_Output") #create directory if it doesn't exist
        os.chdir("Spades_Output") #move to that directory


        #read1 files in input directory 
        fq_files_read1 = []   #list to store forward reads (ending in 1.fq)
        for f in os.listdir(f"../{input_dir}"):
            if f.endswith('1.fq'):
                fq_files_read1.append(f)

        #run SPAdes for each file
        for fq1 in fq_files_read1:     
            fq2 = fq1.replace('1.fq', '2.fq') #find matching reverse read

            fq1_path = f"../{input_dir}/{fq1}"
            fq2_path = f"../{input_dir}/{fq2}"

            #name for output
            base = fq1.replace('1.fq', '')  #create base name output
            outdir = f"{base}"
            os.system(f'mkdir {outdir}') #create directory for each output
        
            #run SPAdes with 2 threads
            os.system(f"spades.py -t 2 --only-assembler -1 {fq1_path} -2 {fq2_path} -o {outdir}") # Hillary's way 
            '''
            os.system("spades.py -t 2 --only-assembler -1 {0} -2 {1} -o Spades_Output".format(fq1_path, fq2_path)) # Josh's way 
            subprocess.run(["spades.py", "-1", fq1_path, "-2", fq2_path, "-o", "Spades_Output"]) # Kimia's way 
            # print(f"finish:{fq1} and {fq2} output in {outdir}")
            '''
        os.chdir('..') #move back to orginal direcotry 

def run_unicycler(): # MAKE CHANGES TO THIS FUNCTION 
    # Unicycler Run 

    #directory to input files
    input_dir = os.path.abspath("artgens") 

    #making output directory 
    if not os.path.isdir("Unicycler_Output"): #make a directory to store SPAdes output if it doesn't already exist 
        os.system("mkdir Unicycler_Output") #create directory if it doesn't exist
    os.chdir("Unicycler_Output") #move to that directory


    #forward read path
    fq1_files =[] #list to store forward reads
    for f in os.listdir(input_dir):
        if f.endswith('1.fq'):
            fq1_files.append(f)

    for fq1 in fq1_files:
        fq2 = fq1.replace('1.fq', '2.fq') #find and match reverse read

        fq1_path_un = os.path.join(input_dir,fq1)
        fq2_path_un = os.path.join(input_dir,fq2)

        base1 = fq1.replace('1.fq', '')  #output folder base name
        outdir1 = f"{base1}"
        os.system(f'mkdir {outdir1}') #create directory for each output
    

        #running unicycler
        os.system(f"unicycler -t 2  -1 {fq1_path_un} -2 {fq2_path_un} -o {outdir1}") #-1 forward read #-2 reverse read  #-o output
        print(f"unicycler finished")


# Kimia's Code - Installing and running QUAST
def install_conda_and_quast():
    current_dir = os.getcwd()  #get thecurrent working directory
    install_dir = os.path.join(current_dir, "quast_install") #the path where miniconda will be installed
    miniconda_path = os.path.join(install_dir, "miniconda3")  #path to miniconda installation
    conda_bin = os.path.join(miniconda_path, "bin", "conda") #path to conda binary
    env_name = "compbio"  #name of conda environment to be created
    python_version = "3.9"

    #If conda is already installed, skip installation
    if not os.path.exists(conda_bin):
        os.makedirs(install_dir, exist_ok=True) #create installation sirectory if not exist
        miniconda_url = "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
        miniconda_script = os.path.join(install_dir, "Miniconda3.sh")
        #download miniconda
        os.system(f"wget -O {miniconda_script} {miniconda_url}")    #-O save file
        #install miniconda
        os.system(f"bash {miniconda_script} -b -p {miniconda_path}") #-b auto accept license  
                                                                       #-P install miniconda to a specific directory
        

    #check if the environment already exists
    env_check = os.popen(f"{conda_bin} env list").read()        #popen runs a shell command and gives back the output

    if env_name in env_check:
        print("conda env already exist")
    else:
        os.system(f"{conda_bin} init bash") #installation conda in bash
        os.system(f"{conda_bin} config --add channels defaults")
        os.system(f"{conda_bin} config --add channels bioconda")
        os.system(f"{conda_bin} config --add channels conda-forge")
        os.system(f"{conda_bin} update conda -y") #update conda
        os.system(f"{conda_bin} create -n {env_name} python={python_version} -y") #create environment
        

    #Install quast
    os.system(f"{conda_bin} run -n {env_name} conda install -c bioconda quast -y")
  
##running quast
#need edit to make sure path directory is working
def run_quast():
    base_dir = os.getcwd()  #current directory
    #define the folder where assemblies are stored  
    spades_dir = os.path.join(base_dir, "spades_output") #spades output folder
    unicycler_dir = os.path.join(base_dir, "unicycler_output")  #unicycler output folder
    modified_ref_dir = os.path.join(base_dir, "Modified_Genomes")    #reference output folder
    quast_output_dir = os.path.join(base_dir, "quast_output")  #output folder for quast


   #get scaffold files
    spades_assembly = os.path.join(spades_dir, "scaffolds.fasta")  #path to spades assembly
    unicycler_assembly = os.path.join(unicycler_dir, "scaffolds.fasta") #path to unicycler assembly

   #finding the reference files
    reference_path = None
    for file in os.listdir(modified_ref_dir):
       if file.endswith(".fna"):  #find modified genome as refrence for quast  #ending in .fna
           reference_path = os.path.join(modified_ref_dir, file)
           break
    #to make sure code can find files    
    if reference_path is None:
        print("reference not found")
        return


    #make output directory 
    os.makedirs(quast_output_dir, exist_ok = True)
    #quast
    quast_run = (
        f"quast.py {spades_dir} {unicycler_dir} "  #spades and unicycler assemblis
        f"-r {reference_path} "      #reference file
        f"-l SPAdes,Unicycler "      #label with spades and unicycler
        f"-o {quast_output_dir}" #output
    )
    print(quast_run)
    os.system(quast_run)
  
# Call the functions in order
if __name__ == "__main__":
    download_genomes()
    generate_and_insert_repeats()
    run_art()
    run_spades()
    run_unicycler()
    install_conda_and_quast()
    run_quast()
