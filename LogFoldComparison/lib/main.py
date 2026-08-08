import sys, os, string
path = os.getcwd()
sys.path.append(os.path.join(path,"lib"))
import tools, foldchange_plot

try:
    import psyco
    psyco.profile()
except:
    pass
    
###############################################################################
# Command line interface
class Interface:
    def __init__(self,options=None):
        self.oValidator = Validator()
        self.cwd = ""
        if __name__ == "__main__":
            self.cwd = ".."
            
        self.options = {
               "-u":"input",        # input folder
               "-o":"output",       # output folder
               "-i":"",             # first pattern
               "-a":"",             # COG column 1
               "-j":"",             # second pattern
               "-b":"",             # COG column 2
               "-t":"",             # COG table
               "-f":"No",           # filter by p-values
                }
        if options:
            self.options.update(options)
            valid = self.oValidator.validate(self.options)
            if valid:
                self.execute()
            else:
                self._main_menu()
        else:
            self._main_menu()

    # Execute selected program
    def execute(self):
        tools.save("\n".join(list(map(lambda item: "%s|%s" % (item[0],str(item[1])),self.options.items()))),os.path.join(self.cwd,"lib","info"))
        self.options['-a'] = int(self.options['-a'])
        self.options['-b'] = int(self.options['-b'])
        if self.options['-f'] != "No":
            self.options['-f'] = float(self.options['-f'])
        title = "%s_vs_%s" % (tools.basename(self.options['-i']),tools.basename(self.options['-j']))
        cog = self._parse_cog()
        transcript1 = self._parse_transcript(self.options['-i'])
        transcript2 = self._parse_transcript(self.options['-j'])
        regulated,pos_inf,neg_inf = self._get_variables(cog,transcript1,transcript2,self.options['-f'])
        svg,output = foldchange_plot.plot([regulated,pos_inf,neg_inf],title,self.options['-f'])
        tools.save(svg,
            os.path.join(self.cwd,self.options['-o'],title+".svg"))
        tools.save(output,
            os.path.join(self.cwd,self.options['-o'],title+".txt"))
        
    def _get_variables(self,cog,transcript1,transcript2,pval_cutoff="No"):
        title1 = tools.basename(self.options['-i'])
        title2 = tools.basename(self.options['-j'])
        regulated = []
        pos_inf = []
        neg_inf = []
        for dataset1,dataset2 in cog:
            dataset1 = list(map(lambda s: s.strip(), dataset1.split("|")))
            dataset2 = list(map(lambda s: s.strip(), dataset2.split("|")))
            try:
                entry1 = list(filter(lambda item: item[0]==dataset1[1], transcript1))[0]
            except:
                continue
            try:
                entry2 = list(filter(lambda item: item[0]==dataset2[1], transcript2))[0]
            except:
                continue
            if not float(entry1[2]) and not float(entry2[2]):
                continue
            if not float(entry1[2]):
                if str(entry2[6])=="NA":
                    continue
                if pval_cutoff != "No" and float(entry2[6]) > pval_cutoff:
                    continue
                neg_inf.append([float(entry2[3]),self._get_product(dataset1[2],dataset1[3]),dataset1[1],entry1[1],[None,float(entry2[6])]]) 
            elif not float(entry2[2]):
                if str(entry1[6])=="NA":
                    continue
                if pval_cutoff != "No" and float(entry1[6]) > pval_cutoff:
                    continue
                pos_inf.append([float(entry1[3]),self._get_product(dataset1[2],dataset1[3]),dataset1[1],entry1[1],[float(entry1[6]),None]])
            else:
                if str(entry1[6])=="NA" or str(entry2[6])=="NA":
                    continue
                if pval_cutoff != "No" and float(entry1[6]) > pval_cutoff and float(entry2[6]) > pval_cutoff:
                    continue
                regulated.append([float(entry1[3]),float(entry2[3]),self._get_product(dataset1[2],dataset1[3]),dataset1[1],entry1[1],
                    [float(entry1[6]),float(entry2[6])]])
        return regulated,pos_inf,neg_inf

    def _get_product(self,gene,product):
        if gene != ".":
            return "%s, %s%s" % (product,gene[0].upper(),gene[1:])
        return product

    def _parse_cog(self):
        def format(ls,w):
            if not ls or not ls[0]:
                return []
            if len(ls) < w:
                ls += ["" for i in range(w-len(ls))]
            elif len(ls) > w:
                ls = ls[:w]
            return ls
        cog = tools.open_text_file(os.path.join(self.cwd,self.options['-u'],self.options['-t']),True,"\t",True)
        width = len(cog[0])
        cog = list(map(lambda item: format(item,width), cog))
        cog = list(filter(lambda item: item, cog))
        cog = list(map(lambda item: [item[self.options['-a']],item[self.options['-b']]], cog))
        return list(filter(lambda item: item[0] and item[0] != "." and item[1] and item[1] != ".", cog))
    
    def _parse_transcript(self,fname):
        transcript = tools.open_text_file(os.path.join(self.cwd,self.options['-u'],fname),True,"\t",True)
        return list(filter(lambda item: len(item) == len(transcript[0]), transcript))[1:]

    # show command prompt interface
    def _main_menu(self):
        response = ''
        while response != "Q":
            print()
            print("Fold change plot of two transcriptoms 2020/05/13")
            print()
            print("Settings for this run:\n")
            print("  I    Transcription file 1\t: " + self.options["-i"])
            print("  A    COG column #1\t\t: " + str(self.options["-a"]))
            print("  J    Transcription file 2\t: " + str(self.options["-j"]))
            print("  B    COG column #2\t\t: " + str(self.options["-b"]))
            print("  T    COG table\t\t: " + self.options["-t"])
            print("  F    P-value cutoff\t\t: " + str(self.options["-f"]))
            print()
            print("Press L-Enter to load the last run options.")
            print("Y to accept these settings, type the letter for one to change or Q to quit")
            print()
            try:
                response = input("?").upper()
                print()
            except:
                continue
            if response == "Q":
                return
            elif response == "Y":
                valid = self.oValidator.validate(self.options)
                if valid:
                    self.execute()
            elif response == "L":
                self._load_options()
            elif response == "I":
                self.options['-i'] = input("Name transcription file 1 in folder '%s'? " % self.options['-u'])
                if not os.path.exists(os.path.join(self.options['-u'],self.options['-i'])):
                    print()
                    print("File %s does not exist!" % os.path.join(self.options['-u'],self.options['-i']))
                    print()
            elif response == "J":
                self.options['-j'] = input("Name transcription file 2 in folder '%s'? " % self.options['-u'])
                if not os.path.exists(os.path.join(self.options['-u'],self.options['-j'])):
                    print()
                    print("File %s does not exist!" % os.path.join(self.options['-u'],self.options['-j']))
                    print()
            elif response == "T":
                self.options['-t'] = input("Name COG table in folder '%s'? " % self.options['-u'])
                if not os.path.exists(os.path.join(self.options['-u'],self.options['-t'])):
                    print()
                    print("File %s does not exist!" % os.path.join(self.options['-u'],self.options['-t']))
                    print()
            elif response == "A":
                try:
                    self.options['-a'] = int(input("Column number in COG table of the 1st dataset? "))
                except:
                    print()
                    print("Column number must be an integer!")
                    print()
            elif response == "B":
                try:
                    self.options['-b'] = int(input("Column number in COG table of the 2nd dataset? "))
                except:
                    print()
                    print("Column number must be an integer!")
                    print()
            elif response == "F":
                v = input("Enter p-value cutoff or 'No'? ")
                try:
                    self.options['-f'] = float(v)
                except:
                    self.options['-f'] = "No"
                
            continue
        
    def _load_options(self):
        path = os.path.join(self.cwd,"lib","info")
        if not os.path.exists(path):
            return
        try:
            options = tools.open_text_file(path,True,"|",True)
        except:
            return
        self.options.update(options)
        
###############################################################################
class Validator:
    def __init__(self):
        self.cwd = ""
        if __name__ == "__main__":
            self.cwd = ".."
        self.prohibited_symbols = [">","<","|",":","\\","/","\"","?","*"]
        
    def validate(self,options,field=""):
        if not field:
            return self.validate_all(options)
        elif field == "-f":
            return True
        elif field in ("-i","-j","-t"):
            if not options[field] or not options['-u']:
                return
            return self.validate_path(os.path.join(self.cwd,options["-u"],options[field]))
        elif field in ("-a","-b"):
            if options[field]=="":
                return
            try:
                options[field]=int(options[field])
            except:
                print()
                print("Column number %s must be a positive integer" % field)
                print()
                return
            if options[field] < 0:
                print()
                print("Column number %s must be a positive integer" % field)
                print()
                return
            cog = tools.open_text_file(os.path.join(self.cwd,options['-u'],options['-t']),True,"\t",True)
            if not cog:
                print()
                print("File %s is corrupted" % os.path.join(self.cwd,options['-u'],options['-t']))
                print()
                return
            if options[field] >= len(cog[0]):
                print()
                print("Column number %s must be a positive integer < %d" % (field,len(cog[0])))
                print()
                return
            return True
        elif field in ("-u","-o"):
            return self.validate_path(os.path.join(self.cwd,options[field]))
        else:
            return
        
    def validate_all(self,options):
        for field in options.keys():
            valid = self.validate(options,field)
            if not valid:
                return
        return True
        
    def validate_path(self,path):
        if not path:
            return
        if os.path.exists(path):
            return path
        return ""
        
    def check_file_name(self,fname):
        if not fname:
            return True
        for symbol in self.prohibited_symbols:
            if fname.find(symbol) > -1:
                return False
        return True
        
###############################################################################

if __name__ == "__main__":
    oInterface = Interface()
