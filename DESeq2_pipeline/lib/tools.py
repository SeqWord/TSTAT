import os, string, math, sys, shelve, time, random
from seq_io import IO
SeqIO = IO()

def save_text_file(path,text):
    f = open(path,"w")
    f.write(str(text))
    f.close()
    return path

def open_text_file(path,inlist=False,separator="",todict=False):
    if not os.path.exists(path):
        return
    try:
        f = open(path)
        data = f.read()
        f.close()
    except:
        return
    if inlist:
        data = list(filter(lambda s: s, data.split("\n")))
        if separator:
            data = list(map(lambda item: list(map(lambda s: s.strip(),item.split(separator))), data))
            if todict:
                if sum(list(map(lambda item: len(item), data))) != 2*len(data):
                    return 
                return dict(zip(list(map(lambda item: item[0], data)),list(map(lambda item: item[1], data))))
    return data

def create_info_list(infolder,extensions=[],datatype="single"):
    if not os.path.exists(infolder):
        return
    if not extensions:
        ls = os.listdir(infolder)
    extensions = list(map(lambda ext: ext.upper(), extensions))
    ls = list(filter(lambda fname: 
        os.path.isfile(os.path.join(infolder,fname)) and fname[fname.rfind("."):].upper() in extensions,
        os.listdir(infolder)))
    if datatype != "single":
        paires = []
        first_files = list(filter(lambda fname: fname.find("_1") != -1, ls))
        ls = list(filter(lambda fname: fname not in first_files, ls))
        second_files = list(map(lambda fname: fname.replace("_1","_2"), first_files))
        for i in range(len(second_files)-1,-1,-1):
            if second_files[i] in ls:
                first_files[i] = "%s,%s" % (first_files[i],second_files[i])
                ls.remove(second_files[i])
            elif datatype=="smart":
                first_files[i] += ", "
        if datatype == "paired-end":
            ls = first_files
        else:
            ls = list(map(lambda fname: fname+", ",ls))+first_files
    return ls
        
def save_info_list(infolder,outfile,extensions=[],datatype="single"):
    ls = create_info_list(infolder,extensions,datatype)
    if not ls:
        return 
    return save_text_file(outfile,"\n".join(ls))
        
def gbk2fasta(gbk_file,outfile,basename=""):
    if not os.path.exists(gbk_file):
        return 
    try:
        gbk = SeqIO.read(gbk_file,"genbank")
    except:
        try:
            gbk = SeqIO.read(open(gbk_file),"genbank")
        except:
            return
    try:
        seq = str(gbk.seq)
        seqname = gbk.description
    except:
        return
    if not basename:
        if seqname.strip():
            basename = seqname
        else:
            basename = "refseq"
    return save_text_file(outfile,">%s\n%s" % (basename,split_seq(seq.upper(),60)))

def gbk2gff(gbk_file,moltypes,outfile,basename):
    gbk = SeqIO.read(gbk_file,"genbank")
    gff = ["##sequence-region %s 1 %d" % (basename,len(gbk.seq))]
    counter = 1
    for i in range(len(gbk.features)):
        if gbk.features[i].type not in moltypes:
            continue
        gene = product = ""
        strand = "+"
        if gbk.features[i].location.strand == -1:
            strand = "-"
        tag = "locus_%d" % counter
        if 'locus_tag' in gbk.features[i].qualifiers:
            tag = gbk.features[i].qualifiers['locus_tag'][0].strip()
        if 'gene' in gbk.features[i].qualifiers:
            gene = gbk.features[i].qualifiers['gene'][0].strip()
        if 'product' in gbk.features[i].qualifiers:
            product = gbk.features[i].qualifiers['product'][0].replace(";",",").strip()
        gff.append("\t".join([basename,"Genbank",gbk.features[i].type,
                              str(int(gbk.features[i].location.start)+1),
                              str(int(gbk.features[i].location.end)),".",strand,"0",
                              "ID=%s-%s;gene=%s;product=%s;locus_tag=%s;gene_id=%d" % (gbk.features[i].type.lower(),tag,gene,product,tag,counter)]
                                       ))
        counter += 1
    return save_text_file(outfile,"\n".join(gff))

def split_seq(seq,n=60):
    start = 0
    stop = n
    length = len(seq)
    if stop >= length:
        return seq
    lines = []
    while stop < length:
        lines.append(seq[start:stop])
        start += n
        stop += n
        if stop >= length:
            lines.append(seq[start:])
            break
    return "\n".join(lines)

def format_number(num,dig,zoom=0):
    return str(int((10**(int(dig)+zoom))*float(num))//float(10**int(dig)))

def scientific_format(num,dig=5):
    if not num:
        return "0.0"
    num = float(num)
    if num >= 1*(10**(-dig)):
        return format_number(num,dig)
    zoom = -math.log(num,10)
    return "%s*10-%d" % (format_number(num,1,zoom),zoom)

def sign(num):
    if num*abs(num) < 0:
        return "-"
    return "+"

def format_string(text,L=30,flg_fill=False):
    if len(text) < L:
        if flg_fill:
            text += (" "*(L-len(text)))
        return text
    else:
        return text[:L-3]+"..."

def dereplicate(ls,sort_fn="",isequal_fnc=None):
    if len(ls) < 2:
        return ls
    if sort_fn:
        ls.sort(eval(sort_fn))
    else:
        ls.sort()
    for i in range(len(ls)-1,0,-1):
        if isequal_fnc:
            if isequal_fnc(ls[i],ls[i-1]):
                del ls[i]
        else:
            if ls[i]==ls[i-1]:
                del ls[i]
    return ls

def validate_path(fname,folder,extensions,markers=[""]):
    if not fname:
        return fname
    if os.path.exists(os.path.join(folder,fname)):
        return os.path.join(folder,fname)
    for ext in extensions:
        for m in markers:
            if os.path.exists(os.path.join(folder,"%s%s%s" % (fname,m,ext))):
                return os.path.join(folder,"%s%s%s" % (fname,m,ext))
    return ""
    
def calculate_pearson_correlation(ls): # ls = [[a1,b1],[a2,b2],...]
    x1 = sum(list(map(lambda item: item[0], ls)))
    x2 = sum(list(map(lambda item: item[0]*item[0], ls)))
    y1 = sum(list(map(lambda item: item[1], ls)))
    y2 = sum(list(map(lambda item: item[1]*item[1], ls)))
    xy = float(sum(list(map(lambda item: item[0]*item[1], ls))))
    d = math.sqrt(len(ls)*x2 - x1**2)*math.sqrt(len(ls)*y2 - y1**2)
    if not d:
        return 0
    return (len(ls)*xy - x1*y1)/d

def calculate_phi_correlation(ctn_tb):
    e = sum(ctn_tb)/4.0
    return abs(ctn_tb[0]-e)/e

def calculate_chi2(ctn_tb):
    return sum(ctn_tb)*float((ctn_tb[0]*ctn_tb[3]-ctn_tb[1]*ctn_tb[2])**2)/(ctn_tb[0]+ctn_tb[1])/(ctn_tb[2]+ctn_tb[3])/(ctn_tb[0]+ctn_tb[2])/(ctn_tb[1]+ctn_tb[3])

def log2(v):
    if v <= 0:
        return 0
    return math.log(abs(v),2)

def log10(p,sign=1.0):
    try:
        p = float(p)
    except:
        return 0
    if not p:
        return 0
    return sign*math.log(p,10)

def average(ls):
    def is_number(v):
        try:
            return float(v)
        except:
            return None
    if not ls:
        return 0
    ls = list(filter(lambda v: v != None, list(map(lambda s: is_number(s), ls))))
    return sum(ls)/len(ls)

goterms = {
        "DNA repair":"DNA replication and repair",
        "DNA replication and repair":"DNA replication and repair",
        "amino aid biosynthesis":"Amino acid and nucleotide base biosynthesis",
        "anabolism":"Anabolism and general biosynthesis",
        "antibiotic":"Antibiotics",
        "biosynthesis":"Anabolism and general biosynthesis",
        "catabolism":"Catabolism and degradation",
        "cell wall biosynthesis":"Cell wall, membrane and transportation",
        "cell wall biosynthesys":"Cell wall, membrane and transportation",
        "degradation":"Catabolism and degradation",
        "enzyme":"Anabolism and general biosynthesis",
        "fatty acid biosynthesis":"Cell wall, membrane and transportation",
        "germination, growth and division":"Cell growth and sporulation",
        "hypothetical":"Hypothetical, phage related and other proteins",
        "membarne":"Cell wall, membrane and transportation",
        "membrane":"Cell wall, membrane and transportation",
        "metabolism":"Anabolism and general biosynthesis",
        "motility, chemotaxis and biofilm formation":"Motility, chemotaxis and biofilm formation",
        "nucleic acid biosynthesis":"Amino acid and nucleotide base biosynthesis",
        "phage associated":"Hypothetical, phage related and other proteins",
        "polyamine biosynthesis":"Stress response and drug resistance",
        "polyamines biosynthesis":"Stress response and drug resistance",
        "polysascharide biosynthesis":"Cell wall, membrane and transportation",
        "resistance":"Stress response and drug resistance",
        "respiration and energy metabolism":"Respiration and energy metabolism",
        "ribosomal":"Protein synthesis",
        "salvage":"Salvage pathways",
        "sporulation":"Cell growth and sporulation",
        "stress response":"Stress response and drug resistance",
        "transcription and protein synthesis":"Protein synthesis",
        "transcriptional factor":"Transcriptional regulation",
        "transcriptional regulator":"Transcriptional regulation",
        "transporter":"Cell wall, membrane and transportation",
        "unidentified":"Hypothetical, phage related and other proteins",
        "vitamin and co-factor biosynthesis":"Vitamin and co-factor biosynthesis"  
    }

goterm_colors = {
        "Amino acid and nucleotide base biosynthesis":"darkturquoise",
        "Anabolism and general biosynthesis":"darkmagenta",
        "Antibiotics":"red",
        "Catabolism and degradation":"darkorange",
        "Cell growth and sporulation":"brown",
        "Cell wall, membrane and transportation":"yellow",
        "DNA replication and repair":"hotpink",
        "Hypothetical, phage related and other proteins":"gainsboro",
        "Motility, chemotaxis and biofilm formation":"aquamarine",
        "Protein synthesis":"lightsalmon",
        "Respiration and energy metabolism":"wheat",
        "Salvage pathways":"peru",
        "Stress response and drug resistance":"green",
        "Transcriptional regulation":"lightskyblue",
        "Vitamin and co-factor biosynthesis":"lightseagreen"
    }

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
        self.oParser = Parser(path)
        fasta = self.oParser.getSequence()
        if flg_ordered_list:
            if len(fasta) > 1:
                seqnames = list(fasta.keys())
                seqnames.sort(key=lambda s: int(s.split("|")[0]))
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
            genes = sorted(genes,key=lambda s: list(map(lambda v: int(v), s.split(".."))))
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
        def get_name(n,gene):
            name = gene['remark'].replace(">","")
            if not name:
                name == "unknown"
            return "%d | %s | [%d..%d]" % (n,name,gene['start'],gene['stop'])
        
        fasta = []
        gene_names = list(genes.keys())
        if len(gene_names) > 1:
            gene_names = sorted(gene_names,key=lambda s: list(map(lambda v: int(v), s.split(".."))))
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
    

##############################################
if __name__ == "__main__":
    pass
    

