import sys, os, string
path = os.getcwd()
sys.path.append(os.path.join(path,"lib"))
import main

###############################################################################
if __name__ == "__main__":
    options = {
               "-u":"input",        # input folder
               "-o":"output",       # output folder
               "-i":"",             # first pattern
               "-a":"",             # COG column 1
               "-j":"",             # second pattern
               "-b":"",             # COG column 2
               "-t":"",             # COG table
               "-f":"No",           # filter by p-values
            }

    arguments = sys.argv[1:]
    if arguments:
        for i in range(0,len(arguments)-1,2):
            key = arguments[i]
            if key not in options:
                raise IOError("Unknown argument " + key + "!")
            if i <= len(arguments)-2:
                options[key] = arguments[i+1]
    
    oMain = main.Interface(options)
