import os, string, math
from functools import reduce
import numpy as np
from scipy import stats

def save(strText,fname=None):
    if not fname:
        return
    ofp = open(fname, "w")
    ofp.write(strText)
    ofp.flush()
    ofp.close()
    return fname

def open_text_file(path,flg_inlist=False,sep="",flg_strip=False):
    if not os.path.exists(path):
        return ""
    f = open(path)
    strText = f.read()
    f.close()
    strText = strText.replace("\"","")
    if flg_inlist:
        strText = strText.split("\n")
        if flg_strip:
            strText = list(map(lambda line: line.strip(), strText))
        if sep:
            strText = list(map(lambda item: item.split(sep), strText))
    return strText

def basename(fname):
    return fname[:fname.rfind(".")]

def matrix(size,fn=float):
    x,y = size
    try:
        x = int(abs(x))
        y = int(abs(y))
    except:
        return []
    return list(list(map(lambda i: y*[fn(0)], range(x))))

def copy_multilevel_ls(ls):
    ls_copy = []
    for item in ls:
        if type(item)==type([]):
            ls_copy.append(copy_multilevel_ls(item))
        else:
            ls_copy.append(item)
    return ls_copy
        
def format_number(num,dig,zoom=0):
    return str(int((10**(dig+zoom))*num)/float(10**dig))

def get_max_value(matrix,members,subset=[]):
    if subset:
        matrix = copy_matrix(matrix,members,subset)
    if len(matrix)<=1:
        return 0
    return max(list(map(lambda v: float(v), reduce(lambda a,b: a+b, matrix))))

def copy_matrix(matrix,members,subset):
    #print "tools:28",len(matrix),len(members),len(subset)
    mx = copy_ls(matrix)
    todelete = list(filter(lambda i: members[i] not in subset, range(len(members))))
    if not todelete:
        return mx
    todelete.reverse()
    for i in todelete:
        del mx[i]
        mx = list(map(lambda row: row[:i]+row[i+1:],mx))
    return mx
    
def copy_ls(ls):
    new_ls = []
    for item in ls:
        if type(item)==type([]):
            new_ls.append(copy_ls(item))
        else:
            try:
                new_ls.append(item.copy())
            except:
                new_ls.append(item)
    return new_ls

def dereplicate(ls,sorting_function=None):
    if len(ls) < 2:
        return ls
    ls.sort(eval(sorting_function))
    for i in range(len(ls)-1,0,-1):
        if ls[i]==ls[i-1]:
            del ls[i]
    return ls

def dereplicate_objects(ls,sorting_function=None):
    if len(ls) < 2:
        return ls
    if sorting_function == None:
        sorting_function = "key=lambda obj: obj.title"
    ls.sort(eval(sorting_function))
    for i in range(len(ls)-1,0,-1):
        if ls[i].title==ls[i-1].title:
            del ls[i]
    return ls
'''
def sort_by_length(a,b):
    return -cmp(len(a),len(b))
'''
def matrix(size,fn=float):
    x,y = size
    try:
        x = int(abs(x))
        y = int(abs(y))
    except:
        return []
    return list(map(lambda i: y*[fn(0)], range(x)))

def calculate_average(ls):
    if not len(ls):
        return 0
    return float(sum(ls))/len(ls)

def calculate_sigma(ls):
    if not len(ls):
        return 0
    n = float(len(ls1))
    x1 = float(sum(ls1))
    x2 = float(sum(list(map(lambda v: v*v, ls1))))
    return x2/n - (x1/n)**2

def calculate_stdev(ls):
    if not len(ls):
        return 0
    return math.sqrt(calculate_sigma(ls))

def contingency_table(data,cutoff=0,moltypes=["CDS","NCRNA"]):
    #9 cell table = [[0,0],[0,0],[0,0],[0,0],...] # [CDS,ncRNA]
    table = []
    table.append(list(map(lambda i: len(list(filter(lambda item: item[0] >= cutoff and item[1] >= cutoff and item[-1].upper()==moltypes[i],data))), 
        range(len(moltypes)))))  # 0: positive coregulation
    table.append(list(map(lambda i: len(list(filter(lambda item: item[0] < -cutoff and item[1] >= cutoff and item[-1].upper()==moltypes[i],data))), 
        range(len(moltypes)))))  # 1: counter-regulation, X negative and Y positive
    table.append(list(map(lambda i: len(list(filter(lambda item: item[0] >= cutoff and abs(item[1]) < cutoff and item[-1].upper()==moltypes[i],data))), 
        range(len(moltypes)))))  # 2: X positively regulated and Y not regulated
    table.append(list(map(lambda i: len(list(filter(lambda item: item[0] <= -cutoff and abs(item[1]) < cutoff and item[-1].upper()==moltypes[i],data))), 
        range(len(moltypes)))))  # 3: X negatively regulated and Y not regulated
    table.append(list(map(lambda i: len(list(filter(lambda item: abs(item[0]) < cutoff and item[1] >= cutoff and item[-1].upper()==moltypes[i],data))), 
        range(len(moltypes)))))  # 4: X not regulated and Y positively regulated
    table.append(list(map(lambda i: len(list(filter(lambda item: abs(item[0]) < cutoff and item[1] <= -cutoff and item[-1].upper()==moltypes[i],data))), 
        range(len(moltypes)))))  # 5: X not regulated and Y negatively regulated
    table.append(list(map(lambda i: len(list(filter(lambda item: abs(item[0]) < cutoff and abs(item[1]) < cutoff and item[-1].upper()==moltypes[i],data))), 
        range(len(moltypes)))))  # 6: X and Y are not regulated
    table.append(list(map(lambda i: len(list(filter(lambda item: item[0] >= cutoff and item[1] < -cutoff and item[-1].upper()==moltypes[i],data))), 
        range(len(moltypes)))))  # 7: counter-regulation, X positive and Y negative
    table.append(list(map(lambda i: len(list(filter(lambda item: item[0] < -cutoff and item[1] < -cutoff and item[-1].upper()==moltypes[i],data))), 
        range(len(moltypes)))))  # 8: negative coregulation
    return table

def calculate_distance(ls): # ls = [[a1,b1],[a2,b2],...]
    d1 = sum(list(map(lambda item: (item[0]-item[1])**2, ls)))
    d2 = sum(list(map(lambda item: 4*max([abs(item[0]),abs(item[1])])**2, ls)))
    if not d2:
        return 0
    return math.sqrt(float(d1)/d2)


def calculate_pearson_correlation(ls): # ls = [[a1,b1],[a2,b2],...]
    query = np.array(list(map(lambda rec: rec[0], ls)))
    sbjct = np.array(list(map(lambda rec: rec[1], ls)))
    if len(query) < 3:
        return 0,1
    pearson_corr,_ = stats.pearsonr(query,sbjct)
    degrees_of_freedom = len(query)-2
    t_stat = pearson_corr * np.sqrt(degrees_of_freedom / (1 - pearson_corr**2))
    p_value = 2 * (1 - stats.t.cdf(np.abs(t_stat), df=degrees_of_freedom))
    return pearson_corr,p_value
    '''
    x1 = sum(list(map(lambda item: item[0], ls)))
    x2 = sum(list(map(lambda item: item[0]*item[0], ls)))
    y1 = sum(list(map(lambda item: item[1], ls)))
    y2 = sum(list(map(lambda item: item[1]*item[1], ls)))
    xy = float(sum(list(map(lambda item: item[0]*item[1], ls))))
    d = math.sqrt(len(ls)*x2 - x1**2)*math.sqrt(len(ls)*y2 - y1**2)
    if not d:
        return 0
    return (len(ls)*xy - x1*y1)/d
    '''

def calculate_phi_correlation(ctn_tb):
    e = sum(ctn_tb)/4.0
    return abs(ctn_tb[0]-e)/e

def calculate_chi2(ctn_tb):
    return sum(ctn_tb)*float((ctn_tb[0]*ctn_tb[3]-ctn_tb[1]*ctn_tb[2])**2)/(ctn_tb[0]+ctn_tb[1])/(ctn_tb[2]+ctn_tb[3])/(ctn_tb[0]+ctn_tb[2])/(ctn_tb[1]+ctn_tb[3])

def calculate_foldchange(ls1,ls2):
    if not ls1 or not ls2:
        return 0,0
    n1 = float(len(ls1))
    n2 = float(len(ls2))
    x1 = float(sum(ls1))
    x2 = float(sum(list(map(lambda v: v*v, ls1))))
    y1 = float(sum(ls2))
    y2 = float(sum(list(map(lambda v: v*v, ls2))))
    avr1 = x1/n1
    avr2 = y1/n2
    sig1 = x2/n1 - (x1/n1)**2
    sig2 = y2/n2 - (y1/n2)**2
    err = ((n1-1)*sig1 + (n2-1)*sig2)*(n1+n2)/(n1+n2-2)/n1/n1
    if not err:
        return 0,0
    t = abs(avr1-avr2)/math.sqrt(err)
    if avr1 <= avr2:
        if not avr1:
            return "+Inf",t
        return avr2/avr1,t
    else:
        if not avr2:
            return "-Inf",t
        return -avr1/avr2,t

# Special methods
def log2fold(v):
    if v == 0:
        return 0
    return math.log(abs(v),2)*int(v+2)/abs(int(v+2))

def log2(v):
    if v <= 0:
        return 0
    return math.log(abs(v),2)

def log10(p,sign=1.0):
    if p <= 0:
        return 0
    return sign*math.log(p,10)

def reverse_logfold(v1,v2,k=2):
    v1 = k**v1
    v2 = k**v2
    return max([v1,v2])/min([v1,v2])

def html2text(s):
    todelete = ['<i>','</i>','<I>','</I>','<sub>','</sub>','<sup>','</sup>','<SUB>','</SUB>',
                '<SUP>','</SUP>','a ', 'an ','A ','An ']
    transtabel = [['&rarr;','&larr;','&harr;','&alpha;','&beta;','&gamma;','&omega;','&zeta;','&pi;','ID\t'],
                  ['=','=','=','alpha','beta','gamma','omega','zeta','p',"'ID\t"]]
    for symbol in todelete:
        s = s.replace(symbol,"")
    for i in range(len(transtabel[0])):
        s = s.replace(transtabel[0][i],transtabel[1][i])
    return s

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
    
pathways = {
        "Biosynthesis":"green",
        "Degradation":"red",
        "Salvage":"yellow",
        "Fatty acid":"hotpink",
        "Gluconeogenesis":"aquamarine",
        "Transport":"darkorange",
    }



