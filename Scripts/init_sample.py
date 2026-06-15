import os

with open(snakemake.input[0]) as files:
    for library in files:
        if not library.startswith("#"):
            tmpLine = library.split('\t') #
            try:
                lib = tmpLine[0]
                fw = tmpLine[1]
                rv = tmpLine[2]
                if lib.lower() == snakemake.wildcards.sample.lower():
                    os.system("Scripts/init_sample.sh "+snakemake.wildcards.PROJECT+" " + lib +" "+fw +" " +rv)
                    exit(0)
            except ValueError:
                print("Error trying to cast: "+ line)
    print("\033[92m There is no entry for Sample: "+ snakemake.wildcards.sample + " in file: " + snakemake.input[0] + " \033[0m")
    print("\033[91m Exiting Meta-Cascabel \033[0m")




