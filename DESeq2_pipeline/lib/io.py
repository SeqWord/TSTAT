import os, sys, math, shelve, string, time, random

###############################################################################
class IO:
    def __init__(self):
        self.oParser = None
        self.valid_fasta_extensions = [".fa",".fasta",".fst",".fna"]
        self.valid_genbank_extensions = [".gbk",".gb"]
        self.valid_seqfiles_extensions = self.valid_fasta_extensions+self.valid_genbank_extensions
    
    def __len__(self):
        if self.oParser:
            return len(self.oParser.getSequence())
        else:
            return 0
    
    #### Collection of save/open functions
    # fasta - FASTA formated collection of sequences    
    def save(self,strText,fname=None):
        if not fname:
            fname = asksaveasfilename([("All files", "*.*")])
        if not fname:
            return
        ofp = open(fname, "w")
        ofp.write(strText)
        ofp.flush()
        ofp.close()
        return fname
    
    def open_text_file(self,path,flg_inlist=False,separator="",flg_strip=False):
        if not os.path.exists(path):
            return ""
        f = open(path)
        strText = f.read()
        f.close()
        if flg_inlist:
            strText = strText.split("\n")
        if flg_inlist and flg_strip:
            strText = [line.strip() for line in strText]
        if flg_inlist and separator:
            strText = list(map(lambda item: item.split(separator), strText))
        return strText
    
    # copy text files
    def copy(self,inpath,outpath):
        if not os.path.exists(inpath):
            return
        try:
            f = open(inpath)
            data = f.read()
            f.close()
            f = open(outpath,"w")
            f.write(data)
            f.close()
        except:
            return
        return outpath
    
    def openGFF(self,path):
        def parse(data):
            data = dict(zip(list(map(lambda item: item.split("=")[0].strip().replace(" ",""),data.split(";"))),
                            list(map(lambda item: item.split("=")[1].strip(),data.split(";")))))
            ls = ['','','']
            if 'locus_tag' in data:
                ls[0] = data['locus_tag']
            if 'gene' in data:
                ls[1] = data['gene']
            if 'product' in data:
                ls[2] = data['product']
            return ls
            
        if not os.path.exists(path):
            return []
        gff = list(filter(lambda item: len(item)==9, self.open_text_file(path,True,"\t",True)))
        data = list(map(lambda item: parse(item[8]), gff))
        data = list(filter(lambda item: item[0], data))
        return dict(zip(list(map(lambda item: item[0], data)),list(map(lambda item: item[1:], data))))
    
    # return dictionary of sequences and the path
    def openFasta(self,path=None,flg_ordered_list=False):
        def sort_seq(a,b):
            return cmp(int(a.split("|")[0]),int(b.split("|")[0]))
        self.oParser = Parser(path)
        fasta = self.oParser.getSequence()
        if flg_ordered_list:
            if len(fasta) > 1:
                seqnames = list(fasta.keys())
                seqnames.sort(sort_seq)
                fasta = [fasta[key] for key in seqnames]
            else:
                fasta = list(fasta.values())
        return fasta,self.oParser.seqlist
    
    # fasta - dictionary {"seqname:seq,...}
    def saveFasta(self,fname,fasta,seqnames=[]):
        if not seqnames:
            seqnames = list(fasta.keys())
        output = []
        for seqname in seqnames:
            output.append(">%s\n%s" % (seqname,fasta[seqname]))
        return self.save("\n".join(output),fname)
    
    # return dictionary of sequences and the path
    def openClustal(self,path=None,filetypes=[]):
        if not filetypes:
            filetypes=[("Alignment files", "*.aln")]        
        self.oParser = Parser(path,filetypes)
        return self.oParser.getSequence(),self.oParser.getPath()
    
    def openGBK(self,path=None,datatype="SEQLIST",filetypes=[]): # datatype - seqlist/dataset/sequence/genemap/fasta/all
        if not filetypes:
            filetypes=[("GBK files", "*.gbk"),
                       ("GBK files", "*.gb")]        
        self.oParser = Parser(path,filetypes)
        if datatype.upper() == "SEQLIST":
            dataset = self.oParser.getDataSet()
            seqname = "%s [%s]" % (dataset['Sequence name'],dataset['Accession'])
            return {seqname:self.oParser.getSequence()},self.oParser.getPath()
        elif datatype.upper() == "DATASET":
            return self.oParser.getDataSet()
        elif datatype.upper() == "SEQUENCE":
            return self.oParser.getSequence()
        elif datatype.upper() == "GENEMAP":
            return self.oParser.getDataSet()['Gene map']
        elif datatype.upper() == "FASTA":
            return ">%s\n%s" % (self.oParser.getName(),self.oParser.getSequence())
        elif datatype.upper() == "ALL":
            return self.oParser.getAll()
        else:
            return
        
    def saveGBK(self,fname,start=0,stop=None,locus_name=""):
        if not self.oParser:
            return
        heading,body,sequence = self.oParser.getGBK_Components(start,stop,locus_name)
        self.save("\n".join([heading,body,self.oParser.format_dna_seq(sequence)]),fname)

    def get_CDS_from_GBK(self,path,ambiguety_threshold=10):
        if not os.path.exists(path):
            return []
        DSet,Seq,path = self.openGBK(path,"ALL")
        oParser = Parser(path)
        return oParser.get_CDC(DSet['Gene map'],Seq)
    
    def openDBFile(self,fname=None, dbkey='$db$', splkey='$suppl$'):
        if not fname or not os.path.exists(fname):
            return
        try:
            f = shelve.open(fname)
            self.oParser = {dbkey:f[dbkey],splkey:{}}
            if splkey in f:
                self.oParser[splkey].update(f[splkey])
            f.close()
            return fname,self.oParser[dbkey],self.oParser[splkey]
        except:
            return None

    def save_genes2fasta(self,gbkpath,outpath,flg_protein=True,lb=0,rb=None,filetypes=None):
        if not filetypes:
            filetypes=[("GBK files", "*.gbk"),
                       ("GBK files", "*.gb")]        
        self.oParser = Parser(gbkpath,filetypes)
        output = self.oParser.genes2fasta(flg_protein,lb,rb)
        if not output:
            return
        self.save(output,outpath)
        return outpath
    
    def getDataSet(self):
        if self.oParser:
            return self.oParser.getDataSet()
        return {}

    def getGeneMap(self):
        if self.oParser:
            return self.oParser.getGeneMap()
        return []
        
    def getSequence(self,seqfile=""):
        if seqfile:
            self.oParser = Parser(seqfile)
        if self.oParser:
            return self.oParser.getSequence()
        return ""
    
    def getName(self):
        if self.oParser:
            return self.oParser.getName()
        return ""
    
    def translate(self,sequences): # seqences = []
        if self.oParser == None:
            self.oParser = Parser()
        sequences = [self.oParser.translate(seq) for seq in sequences]
        self.oParser = None
        return sequences
    
    def clean(self,dirname="",filelist=[]):
        for fname in os.listdir(dirname):
            if fname in filelist or (len(fname)>9 and fname[:-4] in filelist):
                os.remove(os.path.join(dirname,fname))
                continue
            # check for files older than today 
            if (((len(fname)>5 and fname[-4:].upper()==".TMP") or (len(fname)>9 and fname[-8:-3].upper()==".TMP.")) 
                and time.strptime(time.ctime(os.path.getctime(os.path.join(dirname,fname)))).tm_mday != time.localtime().tm_mday):
                os.remove(os.path.join(dirname,fname))
    
    def listdir(self,dirname="",extension=""): # extension without dot
        filelist = []
        for fname in os.listdir(dirname):
            if extension and fname[fname.rfind(".")+1:].upper() != extension.upper():
                continue
            filelist.append(fname)
        return filelist
    
    def random_filename(self,template="",filelist=[]): # template symbols: * for any text symbol; # for any number; ? for any symbol
        i = 1
        while i:
            if len(template) == 0:
                fname = str(random.randint(1,9))
            else:
                fname = ""
                for symbol in template:
                    if symbol == "*":
                        fname += chr(random.randint(97,122))
                    elif symbol == "#":
                        fname += str(random.randint(0,9))
                    elif symbol == "?":
                        code = random.randint(0,1)
                        if code:
                            fname += chr(random.randint(97,122))
                        else:
                            fname += str(random.randint(0,9))
                    elif symbol in (":","<",">","/","\\","|"):
                        fname += "_"
                    else:
                        fname += symbol
            if fname not in filelist:
                return fname
            else:
                i += 1
            if i == 1000:
                print(("System cannot find any vacant file name!"))
                return None

###############################################################################
class Parser:
    def __init__(self, path=None, ftypes=[]):
        # ATTRIBUTES
        self.strSeq = ''
        self.seqlist = []
        self.path = path
        self.DSet = {'Sequence name':'',
                      'Sequence description':'',
                      'Accession':'',
                      'Total sequence length':0,
                      'Locus length':0,
                      'Left border':1,
                      'Frame':0,
                      'Gene map':None,
                      'Path':"",
                      }        
        if self.path and os.path.exists(self.path):
            ext = self.path.split(".")[-1]
            if ext.upper() in ('GBK','GB'):
                self.openGBK()
            elif ext.upper() in ('FA','FAS','FST','FSA','FASTA','FAA','FNN','FNA','ALN'):
                self.openFASTA()
            elif ext.upper() in ('GBF','GBFF'):
                self.openGBFF()
            else:
                pass
        # {'lborder-rborder':[lborder,rborder,dir,gene name,description,remark]}
        self.blast_output = {'Query':'','Sbjct':'','hsps':{}}
        self.codons = {"T":{"T":{"T":"F","C":"F","A":"L","G":"L"},
                       "C":{"T":"S","C":"S","A":"S","G":"S"},
                       "A":{"T":"Y","C":"Y","A":"*","G":"*"},
                       "G":{"T":"C","C":"C","A":"*","G":"W"},
                        },
                   "C":{"T":{"T":"L","C":"L","A":"L","G":"L"},
                       "C":{"T":"P","C":"P","A":"P","G":"P"},
                       "A":{"T":"H","C":"H","A":"Q","G":"Q"},
                       "G":{"T":"R","C":"R","A":"R","G":"R"},
                        },
                   "A":{"T":{"T":"I","C":"I","A":"I","G":"M"},
                       "C":{"T":"T","C":"T","A":"T","G":"T"},
                       "A":{"T":"N","C":"N","A":"K","G":"K"},
                       "G":{"T":"S","C":"S","A":"R","G":"R"},
                        },
                   "G":{"T":{"T":"V","C":"V","A":"V","G":"V"},
                       "C":{"T":"A","C":"A","A":"A","G":"A"},
                       "A":{"T":"D","C":"D","A":"E","G":"E"},
                       "G":{"T":"G","C":"G","A":"G","G":"G"},
                        },
                }

    # METHODS
    
    def __getitem__(self,key):
        if key in self.DSet:
            return self.DSet[key]
        else:
            return ""

    def openFASTA(self,path=""):
        self.seqlist = []
        seqlist = {}
        if not path:
            path = self.path
        objFile = open(path)
        line = objFile.read()
        objFile.close()
        line = line.strip()
        if line[0] != ">":
            self.strSeq = {}
            return
        for symbol in ("\r","\\"):
            line = line.replace(symbol,"")
        line = "\n\r"+line[1:]
        line = line.replace("\n>","\n\r")
        line = line.replace(">","")
        line = line.replace("\n\r",">")
        entries = line.split(">")
        if len(entries) < 2:
            return seqlist
        for i in range(1,len(entries)):
            entry = entries[i]
            data = entry.split("\n")
            if not data[0]:
                continue
            seqlist[data[0]] = ("".join([s.replace(" ","") for s in data[1:]])).upper()
            self.seqlist.append(data[0])
        self.strSeq = seqlist

    # possible modes: 'Get gene map', 'Get sequence', 'Get gene map with sequence'
    def openGBK(self, mode='Get gene map with sequence'):
        file = open(self.path,'r')
        line = file.readline().replace("\r","")
        if line[:5] != 'LOCUS':
            self.openText()
        rb = line.find(" bp")
        line = line[:rb]
        lb = line.rfind(" ")
        self.DSet['Locus length'] = self.DSet['Total sequence length'] = int(line[lb:])
        
        gene = []
        ind = None
        CDS = None
        
        while line:
            #line = file.readline().replace("\r","")
            line = file.readline().strip().replace(">","").replace("<","")
            if "     source          " in line:
                line = line.replace("complement(","")
                try:
                    self.DSet['Left border'] = int(line[string.rfind(line," ")+1:string.rfind(line,"..")])
                except:
                    try:
                        self.DSet['Left border'] = int(line[line.find("join(")+5:line.find("..")])
                    except:
                        self.DSet['Left border'] = 0
            if line.find("DEFINITION  ")==0:
                self.DSet['Sequence name'] = line[12:]
                continue
            if line.find("ACCESSION   ")==0:
                self.DSet['Accession'] = line[12:]
                continue
            #if line[5:9] == 'gene' and mode != 'Get sequence':
            if line[5:8] == 'CDS' and mode != 'Get sequence':
                ind = None
                CDS = 1
                if len(gene) == 6:
                    self.addGene(gene)
                    gene = []
                values = line[21:].split('.')
                if len(values) < 3:
                    continue
                if values[2][0] == ">" or values[2][0] == "<":
                    values[2] = values[2][1:]
                if values[0].find('complement') >= 0:
                    if values[0].find('join') >= 0:
                        try:
                            gene.append(int(values[0][16:]))
                        except:
                            try:
                                gene.append(int(values[0][17:]))
                            except:
                                print((('Error value fot int(): ' + values[0][17:])))
                                return None
                        gene.append(int(values[len(values)-1][:-3]))
                        gene.append('rev')
                    else:   
                        try:
                            gene.append(int(values[0][11:]))
                        except:
                            try:
                                gene.append(int(values[0][12:]))
                            except:
                                print((('Error value fot int(): ' + values[0][12:])))
                                return None
                        gene.append(int(values[2][:-2]))
                        gene.append('rev')
                elif values[0].find('join') >= 0:
                    strand = "dir"
                    try:
                        gene.append(int(values[0][5:]))
                    except:
                        try:
                            gene.append(int(values[0][values[0].find("(")+1:]))
                        except:
                            try:
                                gene.append(int(values[0][values[0].find("complement(")+11:]))
                                strand = "rev"
                            except:
                                print((('Error value fot int(): ' + values[0])))
                                return None
                    try:
                        gene.append(int(values[-1][:values[-1].find(")")]))
                    except:
                        line = file.readline().replace("\r","")
                        while line.find(")") == -1:
                            line = file.readline()
                        gene.append(int(line[line.rfind("..")+2:line.find(")")]))
                    gene.append(strand)
                else:
                    try:
                        gene.append(int(values[0]))
                    except:
                        try:
                            gene.append(int(values[0][1:]))
                        except:
                            print((('Error value for int(): ' + values[0][1:])))
                            return None
                    gene.append(int(values[2]))
                    gene.append('dir')
                for i in range(3):
                    gene.append('')
            elif line[21:22] == r"/" and mode != 'Get sequence' and CDS == 1:
                if line[21:27] == '/gene=' and len(gene) == 6:
                    ind = 3
                    gene[ind] = line[28:-1]
                    if gene[ind] != '' and  gene[ind][-1] == "\"":
                        gene[ind] = gene[ind][:-1]
                        ind = None
                elif line[21:30] == '/product=' and len(gene) == 6:
                    ind = 5
                    gene[ind] = line[31:-1]
                    if gene[ind] != '' and  gene[ind][-1] == "\"":
                        gene[ind] = gene[ind][:-1]
                        ind = None
                elif line[21:27] == '/note=' and len(gene) == 6:
                    ind = 4
                    gene[ind] = line[28:-1]
                    if gene[ind] != '' and gene[ind][-1] == "\"":
                        gene[ind] = gene[ind][:-1]
                        ind = None
                elif line[21:34] == '/translation=':
                    CDS = None
                else:
                    pass
            elif line[:6] == 'ORIGIN' and (mode == 'Get sequence' or mode == 'Get gene map with sequence'):
                ind = None
                if len(gene) == 6:
                    self.addGene(gene)
                self.setSequence(file)
                break
            elif (line == '' or line == '\n') and mode != 'Get sequence' and CDS == 1:
                if len(gene) == 6:
                    self.addGene(gene)
                break
            else:
                if ind and mode != 'Get sequence' and CDS == 1:
                    gene[ind] = gene[ind] + " " + line[21:-1]
                    if gene[ind][-1] == "\"":
                        gene[ind] = gene[ind][:-1]
                        ind = None
                    
        file.close()
        if self.DSet['Total sequence length'] != len(self.strSeq) and len(self.strSeq) != 0:
            self.DSet['Locus length'] = self.DSet['Total sequence length'] = len(self.strSeq)
        return
    
    def openGBFF(self):
        pass

    def addGene(self, gene):
        if self.DSet['Total sequence length'] == 0:
            maxnumlen = 8
        else:
            maxnumlen = len(str(self.DSet['Total sequence length']))
        key = (maxnumlen - len(str(gene[0])))*" " + str(gene[0]) + ".." + str(gene[1])
        if self.DSet['Gene map'] == None:
            self.DSet['Gene map'] = {}
        self.DSet['Gene map'][key] = {}
        subkeys = ('start','stop','direction','name','description','remark')
        for i in range(len(subkeys)):
            self.DSet['Gene map'][key][subkeys[i]] = gene[i]
        return
    
    def e2val(self,e):
        try:
            return float(e)
        except:
            values = string.split(e,"e-")
            if len(values) != 2:
                return None
            if not values[0]:
                values[0] = 1.0
            try:
                return float(values[0]) * (10**(-int(values[1])))
            except:
                return None

    def setSequence(self, file):
        seq = file.read()
        for num in range(10):
            seq = str.replace(seq,str(num),'')
        for symbol in (' ','/','\\','\n'):
            seq = str.replace(seq,symbol,'')
        self.strSeq = str.upper(seq)
        return

    def clear(self):
        self.DSet = {'Sequence name':'',
                        'Accession':'',
                        'Sequence description':'',
                        'Total sequence length':0,
                        'Locus length':0,
                        'Left border':1,
                        'Frame':0,
                        'Gene map':{},
                        }

    # INTERFACE

    # trigger
    def do(self, mode, value=None):
        if mode == 'Set mode':
            self.openGBK(value)
        elif mode == "Import gene map from text file":
            self.openText(value)
        elif mode == "Get line":
            self.showTextViewer(value)
        elif mode == "Set debugger":
            return self.trigger("Set debugger",mode)
        elif mode == "Watch debugger":
            return self.trigger("Watch debugger",mode)
        else:
            print((('Error mode: ' + mode)))

    def getAll(self):
        if self.path:
            return [self.DSet, self.strSeq, self.path]
        else:
            return None

    def getDataSet(self):
        return self.DSet

    def getGeneMap(self):
        if self.path:
            return self.DSet['Gene map']
        else:
            return None
        
    def getSequence(self):
        return self.strSeq
    
    def getPath(self):
        return self.path
    
    def getName(self):
        if self.DSet["Accession"]:
            return self.DSet["Accession"]
        elif self.DSet["Sequence name"]:
            return self.DSet["Sequence name"]
        elif self.DSet["Sequence description"]:
            return self.DSet["Sequence description"]
        else:
            return os.path.basename(self.path)

    def getGBK_Components(self,lb,rb,locus_name="",space=21,width=80,flg_circular=True):
        pre_heading = pre_features = pre_sequence = post_heading = post_features = post_sequence = ''
        
        heading = ["LOCUS       %s" % [s for s in [locus_name,"%s [%d..%d]" % (self.DSet['Accession'],lb,rb)][0] if str(s)]]
        heading[-1] += " "*(40-len(heading[-1])-len(str(rb-lb))) + str(rb-lb) + " bp    DNA     liniar " + self.getTime()
        heading.append("DEFINITION  %s" % self.DSet['Sequence name'])
        heading.append("ACCESSION   %s" % self.DSet['Accession'])
        heading.append("SOURCE      %s" % self.DSet['Sequence name'])
        heading.append("COMMENT     locus start: %d; locus end: %d" % (lb,rb))

        if lb==0:
            lb = 1
        elif lb < 1 and not flg_circular:
            lb = 1
        elif lb < 1 and flg_circulra:
            pre_heading,pre_features,pre_sequence = self.getGBK_Components(len(self.strSeq)+lb,len(self.strSeq),
                    locus_name,space,width)
            lb = 1
        if rb > len(self.strSeq) and not flg_circular:
            rb = len(self.strSeq)
        elif rb > len(self.strSeq) and flg_circular:
            post_heading,post_features,post_sequence = self.getGBK_Components(1,rb-len(self.strSeq)+1,
                    locus_name,space,width)
            rb = len(self.strSeq)

        features = ["FEATURES             Location/Qualifiers"]
        features.append("     source          1.."+str(rb-lb))
        features.append("                     /organism=\""+self.DSet['Sequence name']+"\"")
        if len(features[-1]) > width:
            features[-1] = self.format_string(features[-1],width,space)
                
        genes = list(self.DSet['Gene map'].keys())
        if len(genes) > 1:
            genes.sort()
        for gene in genes:
            try:
                start,stop = [int(s) for s in gene.split("-")]
            except:
                start,stop = [int(s) for s in gene.split("..")]
            if (start < lb and stop <= lb):
                continue
            elif (start >= rb and stop > rb):
                break
            if self.DSet['Gene map'][gene]['direction']=='dir':
                gene_position = "%d..%d" % (start-lb,stop-lb)
            else:
                gene_position = "complement(%d..%d)" % (start-lb,stop-lb)
            features.append("     gene            %s" % gene_position)
            features.append("                     /gene=\"%s\"" % self.DSet['Gene map'][gene]['name'])
            if len(features[-1]) > width:
                features[-1] = self.format_string(features[-1],width,space)
            features.append("                     /db_xref=\"%d..%d\"" % (start,stop))
            if len(features[-1]) > width:
                features[-1] = self.format_string(features[-1],width,space)
            features.append("     CDS             %s" % gene_position) 
            features.append("                     /gene=\"%s\"" % self.DSet['Gene map'][gene]['name'])
            if len(features[-1]) > width:
                features[-1] = self.format_string(features[-1],width,space)
            features.append("                     /product=\"%s\"" % self.DSet['Gene map'][gene]['remark'])
            if len(features[-1]) > width:
                features[-1] = self.format_string(features[-1],width,space)
            aa_seq = self.translate(self.getSequence()[start-1:stop],self.DSet['Gene map'][gene]['direction'])
            features.append("                     /translation=\"%s\"" % self.format_aa_seq(aa_seq))
        return ("\n".join(heading),
                pre_features+"\n".join(features)+post_features,
                pre_sequence+self.getSequence()[lb-1:rb]+post_sequence)
    
    def genes2fasta(self,flg_protein=True,lb=0,rb=None):
        if not self.strSeq:
            return
        genes = list(self.DSet['Gene map'].keys())
        if len(genes) > 1:
            genes.sort(self.sort_genes)
        output = []
        for gene in genes:
            if not gene:
                continue
            try:
                start,stop = [int(s) for s in gene.split("..")]
            except:
                print((gene))
                5/0
            if (start < lb or stop <= lb):
                continue
            elif (rb and (start >= rb or stop > rb)):
                break
            title = "%s, %s (%s) [%d..%d]" % (self.DSet['Gene map'][gene]['name'],
                                            self.DSet['Gene map'][gene]['remark'],
                                            self.DSet['Gene map'][gene]['direction'],
                                            start,stop)
            while title[0] in (" ",","):
                title = title[1:]
            seq = self.substring(start,stop,self.DSet['Gene map'][gene]['direction'])
            if flg_protein:
                seq = self.translate(seq,self.DSet['Gene map'][gene]['direction'])
            output.append(">%s\n%s" % (title,seq))
        return "\n".join(output)
    
    def substring(self,start,stop,strand):
        if strand == "dir":
            if start < 0 and stop >= 0:
                return self.strSeq[start-1:]+self.strSeq[:stop]
            else:
                return self.strSeq[start-1:stop]
        elif strand == "rev":
            if start < 0 and stop >= 0:
                return self.strSeq[start:]+self.strSeq[:stop+1]
            else:
                return self.strSeq[start:stop+1]
        else:
            return
        
    def get_CDC(self,genes,Seq):
        def sort_genes(a,b):
            return cmp(int(a.split("..")[0]),int(b.split("..")[0]))
        
        def get_name(n,gene):
            name = gene['remark'].replace(">","")
            if not name:
                name == "unknown"
            return "%d | %s | [%d..%d]" % (n,name,gene['start'],gene['stop'])
        
        fasta = []
        gene_names = list(genes.keys())
        gene_names.sort(sort_genes)
        for i in range(len(gene_names)):
            gene = genes[gene_names[i]]
            gene_name = get_name(i+1,gene)
            fasta.append(">%s\n%s" % (gene_name,self.translate(Seq[gene['start']-1:gene['stop']],gene['direction'])))
        return fasta
    
    def translate(self,seq,strand='dir'):
        seq = seq.upper()
        start = 0
        aa_seq = ""
        if len(seq) < 3:
            return aa_seq
        if strand == "rev":
            seq = self.reverse_complement(seq)
        codon = seq[start:start+3]
        while start <= len(seq)-3:
            if not codon:
                break
            try:
                aa_seq += self.codons[codon[0]][codon[1]][codon[2]]
            except:
                aa_seq += "X"
            start += 3
            if start >= len(seq)-3:
                break
            codon = seq[start:start+3]
        return aa_seq
    
    def reverse_complement(self,seq):
        seq = seq.upper()
        for s,l in [['A','$'],['T','A'],['$','T'],['C','$'],['G','C'],['$','G']]:
            seq = seq.replace(s,l)
        l = list(seq)
        l.reverse()
        return "".join(l)
    
    def format_aa_seq(self,seq,indend=14,space=21,length=58):
        seq = seq.upper()
        if len(seq) <= length-indend:
            return seq
        i = length-indend
        fseq = [seq[:i]]
        while i < len(seq)-length:
            fseq.append(" "*space + seq[i:i+length])
            i += length
        fseq.append(" "*space + seq[i:])
        return "\n".join(fseq)
    
    def format_dna_seq(self,seq,indend=9,window=60,step=10):
        seq = seq.lower()
        fseq = ["ORIGIN      "]
        i = 1
        while i < len(seq)-window:
            substring = seq[i-1:i+60]
            for j in range(window-step,-1,-step):
                substring = substring[:j]+" "+substring[j:]
            fseq.append(" "*(indend-len(str(i)))+str(i)+substring)
            i += window
        substring = seq[i-1:]
        length = len(substring)
        for j in range(length-length%step,-1,-step):
            substring = substring[:j]+" "+substring[j:]
        fseq.append(" "*(indend-len(str(i)))+str(i)+substring)
        fseq.append("//")
        return "\n".join(fseq)
    
    def format_string(self,seq,width,space):
        j = space
        i = seq.find(" ",j+1)
        border = width
        pos = []
        while i < len(seq):
            if i < 0:
                pos.append(j)
                break
            if i >= border:
                pos.append(j)
                border += width
            j = i
            i = seq.find(" ",j+1)
            
        if pos:
            pos.reverse()
            for p in pos:
                seq = seq[:p] + "\n"+" "*space+seq[p+1:]
        return seq
    
    def getTime(self):
        months = ('JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC')
        year,month,day = time.gmtime()[:3]
        return "%d-%s-%d" % (day,months[month-1],year)
    
    def sort_genes(self,a,b):
        a1,a2 = [int(s) for s in a.split("..")]
        b1,b2 = [int(s) for s in b.split("..")]
        if a1 != b1:
            return a1-b1
        else:
            return a2-b2

##############################################
if __name__ == "__main__":
    pass
    
    
