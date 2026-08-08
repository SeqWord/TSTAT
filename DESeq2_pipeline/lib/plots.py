import os, string, math, random
import tools

###############################################################################
class Plot:
    def __init__(self,experiment,title=""):
        self.experiment = experiment
        self.title = title
        self.oIO = tools.IO()

        # Graphic parameters
        self.width = 800
        self.height = 750
        self.left_margin = 35
        self.top_margin = 50
        self.bottom_margin = 55
        self.right_margin = 100
        self.background = "black"
        self.font_size = 14
        self.plot_width = self.width-2*self.right_margin-2*self.left_margin
        self.plot_height = self.height-2*self.bottom_margin-self.top_margin
    
        self.basic_color = "white"
        self.highlighted_color = "red"
        self.colorscheme = {"protein_coding":[["lightsteelblue","deepskyblue"],["burlywood","coral"]],  # [negative [p1, p2], positive [p1, p2]]
                            "ncrna":[["aquamarine","chartreuse"],["navajowhite","yellow"]]}
        self.positive_color = ["red","orange"]
        self.negative_color = ["chartreuse","deepskyblue"]
        self.basic_dot_size = 2
        self.strock_color = "black"
        self.highlighted_strock_color = "yellow"
        self.strock_width = 1
        
        self.regulated_genes = []
        self.pos_infinity = []
        self.neg_infinity = []

        self.filter_genes = []
        self.highlight_genes = []
        self.outlined_genes = []
        self.source_genome = {}
        
    def _draw_colorscheme(self,X,Y):
        svg = []
        svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
            (X,Y,"white",self.font_size-2,"middle","p"))
        svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
            (X+10,Y,"white",self.font_size-2,"middle","--"))
        svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
            (X+20,Y,"white",self.font_size-2,"middle","0"))
        svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
            (X+30,Y,"white",self.font_size-2,"middle","H"))
        svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
            (X+40,Y,"white",self.font_size-2,"middle","+"))
        svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
            (X+50,Y,"white",self.font_size-2,"middle","p"))
        svg.append("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"%s\" stroke-width=\"1\" />" % 
            (X,Y+10,self.basic_dot_size+2,self._get_color("protein_coding",0,1),"black"))
        svg.append("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"%s\" stroke-width=\"1\" />" % 
            (X+10,Y+10,self.basic_dot_size+2,self._get_color("protein_coding"),"black"))
        svg.append("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"%s\" stroke-width=\"1\" />" % 
            (X+20,Y+10,self.basic_dot_size+2,self.basic_color,"black"))
        svg.append("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"%s\" stroke-width=\"1\" />" % 
            (X+30,Y+10,self.basic_dot_size+2,self.highlighted_color,"black"))
        svg.append("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"%s\" stroke-width=\"1\" />" % 
            (X+40,Y+10,self.basic_dot_size+2,self._get_color("protein_coding",1),"black"))
        svg.append("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"%s\" stroke-width=\"1\" />" % 
            (X+50,Y+10,self.basic_dot_size+2,self._get_color("protein_coding",1,1),"black"))
        svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
            (X+60,Y+13,"white",self.font_size-2,"start","- protein_coding"))
        svg.append("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"%s\" stroke-width=\"1\" />" % 
            (X,Y+25,self.basic_dot_size+2,self._get_color("ncrna",0,1),"black"))
        svg.append("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"%s\" stroke-width=\"1\" />" % 
            (X+10,Y+25,self.basic_dot_size+2,self._get_color("ncrna"),"black"))
        svg.append("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"%s\" stroke-width=\"1\" />" % 
            (X+20,Y+25,self.basic_dot_size+2,self.basic_color,"black"))
        svg.append("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"%s\" stroke-width=\"1\" />" % 
            (X+30,Y+25,self.basic_dot_size+2,self.highlighted_color,"black"))
        svg.append("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"%s\" stroke-width=\"1\" />" % 
            (X+40,Y+25,self.basic_dot_size+2,self._get_color("ncrna",1),"black"))
        svg.append("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"%s\" stroke-width=\"1\" />" % 
            (X+50,Y+25,self.basic_dot_size+2,self._get_color("ncrna",1,1),"black"))
        svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
            (X+60,Y+28,"white",self.font_size-2,"start","- ncrna"))
        return svg

    def _get_color(self,goterm="",i=0,j=0): # i - 0|1 - up or down-regulated; j - 0|1 - p-value higher or lower cutoff
        if not goterm:
            return "deepskyblue"
        if goterm.lower()=="cds":
            goterm = "protein_coding"
        if goterm in self.colorscheme:
            return self.colorscheme[goterm][i][j]
        if goterm in tools.goterms:
            general_term = tools.goterms[goterm]
            return tools.goterm_colors[general_term]
        return "deepskyblue"
    
    def _get_description(self,tag):
        if tag in self.source_genome:
            return tools.format_string(self.source_genome[tag])
        return ""

    def _format_description(self,description,moltype):
        try:
            return tools.format_string("; ".join(list(filter(lambda s: s, [description]+self.source_genome[description]))),80)
        except:
            return "%s; %s" % (moltype,description)

###############################################################################
class VolcanoPlot(Plot):
    def __init__(self,experiment="",title="",filter_file="",highlight_file="",outlined_file="",source_file="",p_5_cutoff=0.05,p_1_cutoff=0.01,inf_cutoff=0.1):
        Plot.__init__(self,experiment,title)
        self.p_1_cutoff = tools.log10(p_1_cutoff,-1.0)
        self.p_5_cutoff = tools.log10(p_5_cutoff,-1.0)
        self.infinity_cutoff = float(inf_cutoff)
        
        self.regulated_genes = []
        self.neg_infinity = []
        self.pos_infinity = []
        self.info = [["p <= %s; fold change >= 2;\nLog2FoldChange\t-Log10(p-value)\tMoltype\tDescription" % tools.scientific_format(p_5_cutoff),
                        ["Up-regulated"],["Down-regulated"]],
                        ["p <= %s; fold change >= 2;\nLog2FoldChange\t-Log10(p-value)\tMoltype\tDescription" % tools.scientific_format(p_1_cutoff),
                        ["Up-regulated"],["Down-regulated"]],
                        ["Positive infinity"],["Negative infinity"],
                    ]
        
        self.occurence_stat = [[0,0],[0,0],[0,0],[0,0]] # [protein_coding,ncrna]
        self.max_fold = 0
        self.min_fold = 0
        self.max_p = 0
        self.max_baseMean = 0
        
        self.x_scale = 0
        self.x_shift = 0
        self.y_scale = 0
        self.y_shift = 0
        self.baseMean_scale = 0
        self.baseMean_shift = 0

        if filter_file and os.path.exists(filter_file) and os.path.isfile(filter_file):
            self.filter_genes = self.oIO.open_text_file(filter_file,True,"\t",True)
            self.filter_genes = list(filter(lambda item: len(item), self.filter_genes))
            self.filter_genes = list(map(lambda item: item[0], self.filter_genes))
        
        if highlight_file and os.path.exists(highlight_file) and os.path.isfile(highlight_file):
            self.highlight_genes = self.oIO.open_text_file(highlight_file,True,"\t",True)
            self.highlight_genes = list(filter(lambda item: len(item), self.highlight_genes))
            self.highlight_genes = list(map(lambda item: item[0], self.highlight_genes))
        
        if outlined_file and os.path.exists(outlined_file) and os.path.isfile(outlined_file):
            self.outlined_genes = self.oIO.open_text_file(outlined_file,True,"\t",True)
            self.outlined_genes = list(filter(lambda item: len(item), self.outlined_genes))
            self.outlined_genes = list(map(lambda item: item[0], self.outlined_genes))
        
        if source_file and os.path.exists(source_file) and os.path.isfile(source_file):
            self.source_genome = self.oIO.openGFF(source_file)
        
    def set(self,data):
        data = list(filter(lambda item: float(item[2]) and item[6] != "NA", data)) # baseMean != 0
        if self.filter_genes:
            data = list(filter(lambda item: item[0] in self.filter_genes, data))
        self.neg_infinity += list(map(lambda item: [item[0],item[1],tools.log10(item[6],-1.0),tools.log2(float(item[2])+1)], 
            list(filter(lambda item: item[9] < self.infinity_cutoff and item[9] < item[8], data)))) # [gene,moltype,p,2*baseMean]
        self.pos_infinity += list(map(lambda item: [item[0],item[1],tools.log10(item[6],-1.0),tools.log2(float(item[2])+1)], 
            list(filter(lambda item: item[8] < self.infinity_cutoff and item[8] < item[9], data)))) # [gene,moltype,p,2*baseMean]
        used_genes = list(map(lambda item: item[0], self.neg_infinity+self.pos_infinity))
        self.regulated_genes += list(map(lambda item: [item[0],item[1],float(item[3]),tools.log10(item[6],-1.0)], # [gene,moltype,fold,p]
            list(filter(lambda item: item[0] not in used_genes, data))))
                    
    def svg(self):    # values = [[[fold1,p1],[fold2,p2]...], [-inf1_baseMean,-inf2_baseMean], [inf1_baseMean,inf2_baseMean]]
        svg = self._svg_template()
        if not self.regulated_genes+self.pos_infinity+self.neg_infinity:
            return svg,""
        svg = svg[:-1]
        self._statistics()
        # Draw line scafold and titles
        svg += self._draw_graphs_and_titles()
        # Draw volcano dots
        svg += self._draw_dots()
        # Draw occurence statistics
        svg += self._draw_occurence_stat()
        # Draw negative infinity 
        svg += self._draw_infinity_column(self.neg_infinity,.5,0,5) 
        # Draw positive infinity 
        svg += self._draw_infinity_column(self.pos_infinity,1.5,1,-5) 
        # Draw Legends
        svg += self._draw_legend()
        
        if self.experiment:
            svg.append("</g>")
        svg.append("</svg>")
        
        return [svg, ("%s\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s" % 
                    (self.info[0][0],"\n".join(self.info[0][1]),
                    "\n".join(self.info[0][2]),self.info[1][0],
                    "\n".join(self.info[1][1]),
                    "\n".join(self.info[1][2]),
                    "\n".join(self.info[2]),
                    "\n".join(self.info[3])
                    ))]

    def _draw_graphs_and_titles(self):
        svg = []
        if self.title:
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                (self.left_margin+self.plot_width/2,self.height-5,"orange",self.font_size+4,"middle",self.title))
        if self.experiment:
            svg.insert(1,"<g id=\"%s\">" % self.experiment)
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                (3*self.left_margin,self.top_margin/2,"white",self.font_size+2,"start",self.experiment))

        # X-scale
        for i in range(int(self.max_fold-self.min_fold)):
            x = 2*self.left_margin+(self.x_shift+i)*self.x_scale
            svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
                (x,self.top_margin+self.plot_height,x,self.top_margin+self.plot_height+5,"white",1))
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%d</text>\"" % 
                (x,self.top_margin+self.plot_height+7+self.font_size,"white",self.font_size,"middle",int(self.min_fold)+i))
        # vertical line -1
        i = -1-int(self.min_fold)
        x = 2*self.left_margin+(self.x_shift+i)*self.x_scale
        svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
            (x,self.top_margin,x,self.top_margin+self.plot_height+5,"grey",1))
        # vertical line 1
        i = 1-int(self.min_fold)
        x = 2*self.left_margin+(self.x_shift+i)*self.x_scale
        svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
            (x,self.top_margin,x,self.top_margin+self.plot_height+5,"grey",1))
        
        # Y-scale
        y_mark_spacing = int(self.max_p/10)
        for i in range(int(self.max_p)):
            y = self.top_margin+(self.y_shift+i)*self.y_scale
            svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
                (2*self.left_margin-5,y,2*self.left_margin,y,"white",1))
            if not y_mark_spacing or (int(self.max_p)-i)%y_mark_spacing == 0:
                svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%d</text>\"" % 
                    (2*self.left_margin-10,y,"white",self.font_size,"end",int(self.max_p)-i))
        # p 0.05 horizontal line
        y = self.top_margin+(0.5+self.max_p-self.p_5_cutoff)*self.y_scale
        svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
            (2*self.left_margin,y,2*self.left_margin+self.plot_width,y,"grey",1))
        # p 0.01 horizontal line
        y = self.top_margin+(0.5+self.max_p-self.p_1_cutoff)*self.y_scale
        svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
            (2*self.left_margin,y,2*self.left_margin+self.plot_width,y,"white",1))
            
        # baseMean-scale
        z_mark_spacing = int(self.max_baseMean/10)
        if self.max_baseMean:
            for i in range(int(self.max_baseMean)):
                y = self.top_margin+(self.baseMean_shift+i)*self.baseMean_scale
                svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
                    (2*self.left_margin+self.plot_width+self.right_margin-5,y,2*self.left_margin+self.plot_width+self.right_margin+5,y,"white",1))
                if not z_mark_spacing or (int(self.max_baseMean)-i)%z_mark_spacing == 0:
                    svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%d</text>\"" % 
                        (2*self.left_margin+self.plot_width+self.right_margin+10,y,"white",self.font_size,"start",int(self.max_baseMean)-i))
        return svg
    
    def _draw_dots(self):
        svg = []
        # Add volcano dots
        for i in range(len(self.regulated_genes)):
            gene,moltype,fold,p = self.regulated_genes[i]
            description = self._format_description(gene,moltype)
            x = 2*self.left_margin+(fold-self.min_fold+0.5)*self.x_scale
            y = self.top_margin+(.5+self.max_p-p)*self.y_scale
            r = self.basic_dot_size
            color = self.basic_color
            reg = dm = 0
            if p >= self.p_5_cutoff and abs(fold) >= 1:
                r += 2
                if fold < 0:
                    self.info[0][2].append("%f\t%f\t%s\t%s" % (fold,p,moltype,description))
                    if moltype.lower() in ("protein_coding","cds"):
                        self.occurence_stat[1][0] += 1
                    else:
                        self.occurence_stat[1][1] += 1
                else:
                    reg = 1
                    self.info[0][1].append("%f\t%f\t%s\t%s" % (fold,p,moltype,description))
                    if moltype.lower() in ("protein_coding","cds"):
                        self.occurence_stat[2][0] += 1
                    else:
                        self.occurence_stat[2][1] += 1
                color = self._get_color(moltype,reg,dm)
            if p >= self.p_1_cutoff and abs(fold) >= 1:
                r += 2
                dm = 1
                if fold < 0:
                    self.info[1][2].append("%f\t%f\t%s\t%s" % (fold,p,moltype,description))
                    if moltype.lower() in ("protein_coding","cds"):
                        self.occurence_stat[1][0] -= 1
                        self.occurence_stat[0][0] += 1
                    else:
                        self.occurence_stat[1][1] -= 1
                        self.occurence_stat[0][1] += 1
                else:
                    self.info[1][1].append("%f\t%f\t%s\t%s" % (fold,p,moltype,description))
                    if moltype.lower() in ("protein_coding","cds"):
                        self.occurence_stat[2][0] -= 1
                        self.occurence_stat[3][0] += 1
                    else:
                        self.occurence_stat[2][1] -= 1
                        self.occurence_stat[3][1] += 1
                color = self._get_color(moltype,reg,dm)
            if gene in self.highlight_genes:
                color = self.highlighted_color 
            strock_color = self.strock_color
            strock_width = self.strock_width
            if gene in self.outlined_genes:
                strock_color = self.highlighted_strock_color
                strock_width += 1
            svg.append(("<a xlink:title=\"%s\">" % description)+"<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"%s\" stroke-width=\"%d\" /></a>" % 
                        (x,y,r,color,strock_color,strock_width))
        return svg
    
    def _draw_occurence_stat(self):
        svg = []
        # Stat values
        y1 = self.top_margin+(0.5+self.max_p-self.p_1_cutoff)*self.y_scale-5
        y2 = self.top_margin+(0.5+self.max_p-self.p_5_cutoff)*self.y_scale-5
        if self.occurence_stat[0]:
            if self.occurence_stat[0][0] and self.occurence_stat[0][1]:
                svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                    (2*self.left_margin+5,y1,"yellow",self.font_size-2,"start",
                    "%dc+%dnc" % (self.occurence_stat[0][0],self.occurence_stat[0][1])))
            elif self.occurence_stat[0][0]:
                svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                    (2*self.left_margin+5,y1,"yellow",self.font_size-2,"start",
                    "%d cds" % (self.occurence_stat[0][0])))
            else:
                svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                    (2*self.left_margin+5,y1,"yellow",self.font_size-2,"start",
                    "%d ncrna" % self.occurence_stat[0][1]))
        if self.occurence_stat[1]:
            if self.occurence_stat[0][0]+self.occurence_stat[1][0] and self.occurence_stat[0][1]+self.occurence_stat[1][1]:
                svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                    (2*self.left_margin+5,y2,"yellow",self.font_size-2,"start",
                    "%dc+%dnc" % (self.occurence_stat[0][0]+self.occurence_stat[1][0],self.occurence_stat[0][1]+self.occurence_stat[1][1])))
            elif self.occurence_stat[0][0]+self.occurence_stat[1][0]:
                svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                    (2*self.left_margin+5,y2,"yellow",self.font_size-2,"start",
                    "%d cds" % (self.occurence_stat[0][0]+self.occurence_stat[1][0])))
            else:
                svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                    (2*self.left_margin+5,y2,"yellow",self.font_size-2,"start",
                    "%d ncrna" % (self.occurence_stat[0][1]+self.occurence_stat[1][1])))
        
        if self.occurence_stat[2]:
            if self.occurence_stat[2][0]+self.occurence_stat[3][0] and self.occurence_stat[2][1]+self.occurence_stat[3][1]:
                svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                    (2*self.left_margin+self.plot_width-5,y2,"yellow",self.font_size-2,"end",
                    "%dc+%dnc" % (self.occurence_stat[2][0]+self.occurence_stat[3][0],self.occurence_stat[2][1]+self.occurence_stat[3][1])))
            elif self.occurence_stat[2][0]+self.occurence_stat[3][0]:
                svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                    (2*self.left_margin+self.plot_width-5,y2,"yellow",self.font_size-2,"end",
                    "%d cds" % (self.occurence_stat[2][0]+self.occurence_stat[3][0])))
            else:
                svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                    (2*self.left_margin+self.plot_width-5,y2,"yellow",self.font_size-2,"end",
                    "%d ncrna" % (self.occurence_stat[2][1]+self.occurence_stat[3][1])))
        if self.occurence_stat[3]:
            if self.occurence_stat[3][0] and self.occurence_stat[3][1]:
                svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                    (2*self.left_margin+self.plot_width-5,y1,"yellow",self.font_size-2,"end",
                    "%dc+%dnc" % (self.occurence_stat[3][0],self.occurence_stat[3][1])))
            elif self.occurence_stat[3][0]:
                svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                    (2*self.left_margin+self.plot_width-5,y1,"yellow",self.font_size-2,"end",
                    "%d cds" % (self.occurence_stat[3][0])))
            else:
                svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                    (2*self.left_margin+self.plot_width-5,y1,"yellow",self.font_size-2,"end",
                    "%d ncrna" % (self.occurence_stat[3][1])))
        return svg
    
    def _draw_infinity_column(self,infinity_ls,shift,reg,dent):  # shift - 0.5|1.5; reg - 0|1
        sign = [-1.0,1.0]
        svg = []
        for i in range(len(infinity_ls)):
            gene,moltype,p,v = infinity_ls[i]
            description = self._format_description(gene,moltype)
            x = 2*self.left_margin+self.plot_width+shift*self.right_margin+random.randint(-10*self.basic_dot_size,10*self.basic_dot_size)
            y = self.top_margin+(.5+self.max_baseMean-v)*self.baseMean_scale
            r = self.basic_dot_size
            color = self.basic_color
            dm = 0
            if p >= self.p_5_cutoff:
                dm = 1
            if v >= 2:
                r += 2
                if v >= 4:
                    r += 2
                color = self._get_color(moltype,reg,dm)
                self.info[3-reg].append("\t".join([str(v*sign[reg]),str(p),moltype,description]))
            if gene in self.highlight_genes:
                color = self.highlighted_color 
            strock_color = self.strock_color
            strock_width = self.strock_width
            if gene in self.outlined_genes:
                strock_color = self.highlighted_strock_color
                strock_width += 1
            svg.append(("<a xlink:title=\"%s\">" % description)+"<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"%s\" stroke-width=\"%d\" /></a>" % 
                        (x,y,r,color,strock_color,strock_width))
        cds_count = len(list(filter(lambda item: item[1].lower() in ("protein_coding","cds") and item[3] >= 2, infinity_ls)))
        ncrna_count = len(list(filter(lambda item: item[1].lower() in ("ncrna","misc_rna","rna") and item[3] >= 2, infinity_ls)))
        if cds_count or ncrna_count:
            if cds_count and ncrna_count:
                legend = "%dc+%dnc" % (cds_count,ncrna_count)
            elif cds_count:
                legend = "%d cds" % cds_count
            else:
                legend = "%d ncrna" % ncrna_count
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                (5+2*self.left_margin+self.plot_width+(shift-dent/10.0)*self.right_margin/2+dent-dent,self.top_margin+self.plot_height+15,"yellow",self.font_size-2,"start",legend))
        return svg
    
    def _draw_legend(self):
        svg = []
        svg += self._draw_colorscheme(2*self.left_margin+self.plot_width+self.right_margin/2,self.top_margin+self.plot_height+70)
        return svg
    
    def _svg_template(self):
        svg = ["<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" viewbox=\"0 0 %d %d\">" % (self.width,self.height)]
        # background
        svg.append("<rect x=\"0\" y=\"0\" fill=\"%s\" width=\"%d\" height=\"%d\"/>" % (self.background,self.width,self.height))
        # left main vertical axis
        svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
            (2*self.left_margin,self.top_margin,2*self.left_margin,self.height-self.bottom_margin,"white",2))
        # right main vertical axis
        svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
            (self.width-2*self.right_margin,self.top_margin,self.width-2*self.right_margin,self.height-self.bottom_margin,"white",2))
        # right baseMean axis
        svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
            (self.width-self.right_margin,self.top_margin,self.width-self.right_margin,self.height-self.bottom_margin,"white",1))
        # horizontal axis
        svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
            (self.left_margin,self.height-2*self.bottom_margin,self.width-10,self.height-2*self.bottom_margin,"white",2))
        # vertical axis titles
        svg.append(("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" "+
            "transform=\"rotate(%f %d,%d)\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>") %
            (self.left_margin,self.top_margin+(self.height-2*self.bottom_margin-self.top_margin)/2,-90,self.left_margin,self.top_margin+(self.height-2*self.bottom_margin-self.top_margin)/2,
            "white",self.font_size,"middle","-Log10(p-value)"))
        svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
            (self.width-self.right_margin,self.top_margin-1.5*self.font_size,"white",self.font_size,"middle","Log2(baseMean+1)"))
        # horizontal axis titles
        svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
            (self.left_margin+(self.width-2*self.right_margin-2*self.left_margin)/2,self.height-self.bottom_margin,"white",self.font_size,"middle","Log2FoldChange"))
        svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
            (self.width-1.5*self.right_margin,self.height-self.bottom_margin,"white",self.font_size,"middle","-Inf"))
        svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
            (self.width-.5*self.right_margin,self.height-self.bottom_margin,"white",self.font_size,"middle","Inf"))

        svg.append("</svg>")
        return svg
    
    def _statistics(self):
        self.max_fold = max(list(map(lambda ls: float(ls[2]), self.regulated_genes)))
        self.min_fold = min(list(map(lambda ls: float(ls[2]), self.regulated_genes)))
        self.max_p = max(list(map(lambda ls: float(ls[3]), self.regulated_genes)))
        
        self.x_scale = float(self.plot_width)/(self.max_fold-self.min_fold+1.0)
        self.x_shift = 0.5+abs(self.min_fold)-abs(int(self.min_fold))
        self.y_scale = float(self.plot_height)/(self.max_p+0.5)
        self.y_shift = 0.5+self.max_p-int(self.max_p)

        self.max_baseMean = self.baseMean_scale = self.baseMean_shift = 0
        if self.pos_infinity or self.neg_infinity:
            self.max_baseMean = max(list(map(lambda ls: float(ls[3]), self.pos_infinity+self.neg_infinity)))
            self.baseMean_scale = float(self.plot_height)/(self.max_baseMean+0.5)
            self.baseMean_shift = 0.5+abs(self.max_baseMean)-abs(int(self.max_baseMean))

###############################################################################
class ExpressionPlot(Plot):
    def __init__(self,strain1,strain2,title="",filter_file="",highlight_file="",outlined_file="",source_file="",p_cutoff=0.05,infinity_cutoff=0.001):
        self.strain1 = strain1
        self.strain2 = strain2
        self.generic_title = "Log2(baseMean+1)"
        Plot.__init__(self,"%s_vs_%s" % (self.strain1,self.strain2),title)
        
        if filter_file and os.path.exists(filter_file) and os.path.isfile(filter_file):
            self.filter_genes = self.oIO.open_text_file(filter_file,True,"\t",True)
            self.filter_genes = list(filter(lambda item: len(item), self.filter_genes))
            self.filter_genes = list(map(lambda item: item[0], self.filter_genes))
        
        if highlight_file and os.path.exists(highlight_file) and os.path.isfile(highlight_file):
            self.highlight_genes = self.oIO.open_text_file(highlight_file,True,"\t",True)
            self.highlight_genes = list(filter(lambda item: len(item), self.highlight_genes))
            self.highlight_genes = list(map(lambda item: item[0], self.highlight_genes))
        
        if outlined_file and os.path.exists(outlined_file) and os.path.isfile(outlined_file):
            self.outlined_genes = self.oIO.open_text_file(outlined_file,True,"\t",True)
            self.outlined_genes = list(filter(lambda item: len(item), self.outlined_genes))
            self.outlined_genes = list(map(lambda item: item[0], self.outlined_genes))
        
        if source_file and os.path.exists(source_file) and os.path.isfile(source_file):
            self.source_genome = self.oIO.openGFF(source_file)
        
        self.p_cutoff = p_cutoff
        self.infinity_cutoff = float(infinity_cutoff)
        
        self.categories = ["Up-regulated in %s compared to %s" % (self.strain1,self.strain2),
                    "Down-regulated in %s compared to %s" % (self.strain1,self.strain2),
                    "Expressed in %s" % self.strain1,"Expressed in %s" % self.strain2]
        self.info = dict(zip(self.categories,list(map(lambda i: [], range(len(self.categories))))))
        
        self.maxX = 0
        self.minX = 0
        self.maxY = 0
        self.minY = 0
        self.maxZ = 0
        self.minZ = 0
        self.x_scale = 0
        self.x_shift = 0
        self.y_scale = 0
        self.y_shift = 0
        self.side_scale = 0
        
        self.x_avr = 0
        self.y_avr = 0
        self.expression_degree = 0
        self.pearson = 0
 
        self.regulated_genes = []   # [gene,moltype,p,baseMean1,baseMean2]
        self.pos_infinity = []      # [gene,moltype,baseMean2]
        self.neg_infinity = []      # [gene,moltype,baseMean1]

    def set(self,data):
        data = list(filter(lambda item: float(item[2]) and item[6] != "NA", data)) # baseMean != 0
        if self.filter_genes:
            data = list(filter(lambda item: item[0] in self.filter_genes, data))
        self.neg_infinity += list(map(lambda item: [item[0],item[1],float(item[6]),tools.log2(2.0*float(item[8])+1)], 
            list(filter(lambda item: item[9] < self.infinity_cutoff and item[9] < item[8], data)))) # [gene,moltype,p,2*baseMean]
        self.pos_infinity += list(map(lambda item: [item[0],item[1],float(item[6]),tools.log2(2.0*float(item[9])+1)], 
            list(filter(lambda item: item[8] < self.infinity_cutoff and item[8] < item[9], data)))) # [gene,moltype,p,2*baseMean]
        used_genes = list(map(lambda item: item[0], self.neg_infinity+self.pos_infinity))
        self.regulated_genes += list(map(lambda item: [item[0],item[1],2**abs(float(item[3])),float(item[6]),
            tools.log2(float(item[8])+1),tools.log2(float(item[9])+1)], # [gene,moltype,foldchange,p,baseMean1,baseMean2]
            list(filter(lambda item: item[0] not in used_genes, data))))
                    
    def svg(self): 
        svg = self._svg_template(self.strain1,self.strain2)
        if not self.regulated_genes+self.pos_infinity+self.neg_infinity:
            return svg,""
        svg = svg[:-1]
        self._statistics()
        # Draw line scafold and titles
        svg += self._draw_graphs_and_titles()
        # Draw volcano dots
        svg += self._draw_dots()
        # Draw trend line
        svg += self._draw_trend()
        # Draw occurence statistics
        svg += self._draw_occurence_stat()
        # Draw negative infinity 
        svg += self._draw_infinity_column(self.neg_infinity,.5,0,5) 
        # Draw positive infinity 
        svg += self._draw_infinity_column(self.pos_infinity,1.5,1,-5)
        # Draw Legends
        svg += self._draw_legend()
        
        if self.experiment:
            svg.insert(1,"<g id=\"%s\">" % self.experiment)
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s<tspan font-size=\"%f\">%s</tspan></text>\"" % 
                (3*self.left_margin,self.top_margin/2,"white",self.font_size+2,"start",self.experiment,self.font_size,
                " Pearson Corr. = %s; expr. degree = %s" % (tools.format_number(self.pearson,3),tools.format_number(self.expression_degree,1))))
            svg.append("</g>")
        svg.append("</svg>")
        output = []
        for i in range(2):
            key_cat = self.categories[i]
            self.info[key_cat] = sorted(list(self.info[key_cat]),key=lambda ls: ls[-3]**2+ls[-2]**2,reverse=True)
            output.append("%s\n\t%s" % (key_cat,"\n\t".join(list(map(lambda item: "\t".join(list(map(lambda v: str(v), item))), self.info[key_cat])))))
        for i in range(2,4,1):
            key_cat = self.categories[i]
            self.info[key_cat] = sorted(list(self.info[key_cat]), key=lambda ls: ls[-1],reverse=True)
            output.append("%s\n\t\t%s" % (key_cat,"\n\t\t".join(list(map(lambda item: "\t".join(list(map(lambda v: str(v), item))), self.info[key_cat])))))
        output.insert(0,"\n".join(["# Locus tag","# MolType","# Fold change","# p-value","# Log2(first normalized count)","# Log2(second normalized count)"]))
        return svg,"\n".join(output)

    def _draw_graphs_and_titles(self):
        svg = []
        if self.title:
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                (self.left_margin+self.plot_width/2,self.height-5,"orange",self.font_size+4,"middle",self.title))
        # X-scale
        for i in range(int(self.maxX-self.minX)+1):
            x = 2*self.left_margin+(self.x_shift+i)*self.x_scale
            if x >= 2*self.left_margin+self.plot_width:
                break
            svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
                (x,self.top_margin+self.plot_height,x,self.top_margin+self.plot_height+5,"white",1))
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%d</text>\"" % 
                (x,self.top_margin+self.plot_height+7+self.font_size,"white",self.font_size,"middle",int(self.minX)+i))
        # vertical line 0
        i = -int(self.minX)
        x = 2*self.left_margin+(self.x_avr-self.minX+0.5)*self.x_scale
        svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
            (x,self.top_margin,x,self.top_margin+self.plot_height+5,"grey",1))
        
        # Y-scale
        for i in range(int(self.maxY-self.minY)+1):
            y = self.top_margin+(.5+self.maxY-int(self.maxY)+i)*self.y_scale
            if y >= self.top_margin + self.plot_height:
                break
            svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
                (2*self.left_margin-5,y,2*self.left_margin,y,"white",1))
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%d</text>\"" % 
                (2*self.left_margin-10,y,"white",self.font_size,"end",int(self.maxY)-i))
        # horizontal line 0
        y = self.top_margin+(.5+self.maxY-self.y_avr)*self.y_scale
        svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
            (2*self.left_margin,y,2*self.left_margin+self.plot_width,y,"grey",1))
        # Side-scale
        if self.maxZ+self.minZ:
            side_axis_shift = 0.5+abs(self.minZ)-abs(int(self.minZ))
            if self.maxZ > 1 and self.minZ < -1:
                # side horizontal line 0
                y = self.top_margin+(self.maxZ+.5)*self.side_scale
                svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
                    (2*self.left_margin+self.plot_width,y,2*self.left_margin+self.plot_width+2*self.right_margin,y,"grey",1))
            for i in range(int(self.maxZ-self.minZ)+1):
                y = self.top_margin+(.5+self.maxZ-int(self.maxZ)+i)*self.side_scale
                if y >= self.top_margin + self.plot_height:
                    break
                svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
                    (2*self.left_margin+self.plot_width+self.right_margin-5,y,2*self.left_margin+self.plot_width+self.right_margin+5,y,"white",1))
                svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%d</text>\"" % 
                    (2*self.left_margin+self.plot_width+self.right_margin+10,y,"white",self.font_size,"start",int(self.maxZ)-i))
        return svg
    
    def _draw_dots(self):
        svg = []
        # Add dots
        for i in range(len(self.regulated_genes)):
            gene,moltype,foldchange,p,RPKM1,RPKM2 = self.regulated_genes[i]
            description = self._format_description(gene,moltype)
            x = 2*self.left_margin+(RPKM1-self.minX+0.5)*self.x_scale
            y = self.top_margin+(.5+self.maxY-RPKM2)*self.y_scale
            r = self.basic_dot_size
            color = self.basic_color
            if foldchange >= 2:
                reg = dm = 0
                if RPKM2 > RPKM1:
                    reg = 1
                if p <= self.p_cutoff:
                    dm = 1
                color = self._get_color(moltype,reg,dm)
                if RPKM1 < RPKM2:
                    self.info[self.categories[0]].append(self.regulated_genes[i]+[tools.format_string(description,100)])
                else:
                    self.info[self.categories[1]].append(self.regulated_genes[i]+[tools.format_string(description,100)])
                if foldchange >= 4:
                    r += 2
                    if foldchange >= 8:
                        r += 2
            if gene in self.highlight_genes:
                color = self.highlighted_color 
            strock_color = self.strock_color
            strock_width = self.strock_width
            if gene in self.outlined_genes:
                strock_color = self.highlighted_strock_color
                strock_width += 1
            svg.append(("<a xlink:title=\"%s\">" % description)+("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"%s\" stroke-width=\"%d\" /></a>" % 
                        (x,y,r,color,strock_color,strock_width)))
        return svg
    
    def _draw_trend(self):
        svg = []
        # trend line
        if self.expression_degree:
            if self.y_avr:
                slop = self.x_avr/self.y_avr
            elif self.y_avr - self.minY:
                slop = (self.x_avr-self.minX)/(self.y_avr-self.minY)
            else:
                return svg
            if self.expression_degree < 0:
                x = 2*self.left_margin+(self.maxX-self.minX+1.0)*self.x_scale
                y = self.top_margin+(self.maxY-self.maxX/slop)*self.y_scale
            else:
                x = 2*self.left_margin+(self.maxY*slop-self.minX+1.0)*self.x_scale
                y = self.top_margin
            svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
                (2*self.left_margin,self.top_margin+self.plot_height,x,y,"red",1))
        return svg
    
    def _draw_occurence_stat(self):
        svg = []
        # Stat values
        if self.info[self.categories[1]]:
            cds_n = len(list(filter(lambda ls: ls[1].lower() in ("protein_coding","cds"), self.info[self.categories[0]])))
            rna_n = len(list(filter(lambda ls: ls[1].lower() in ("ncrna","misc_rna","rna"), self.info[self.categories[0]])))
            if rna_n:
                svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                    (2*self.left_margin+self.plot_width-15,self.top_margin+self.plot_height-5,"yellow",self.font_size-2,"end","%dc+%dnc" % (cds_n,rna_n)))
            else:
                svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                    (2*self.left_margin+self.plot_width-15,self.top_margin+self.plot_height-5,"yellow",self.font_size-2,"end","%d CDS" % cds_n))
        
        if self.info[self.categories[0]]:
            cds_n = len(list(filter(lambda ls: ls[1].lower() in ("protein_coding","cds"), self.info[self.categories[1]])))
            rna_n = len(list(filter(lambda ls: ls[1].lower() in ("ncrna","misc_rna","rna"), self.info[self.categories[1]])))
            if rna_n:
                svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                    (2*self.left_margin+15,self.top_margin+15,"yellow",self.font_size-2,"start","%dc+%dnc" % (cds_n,rna_n)))
            else:
                svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                    (2*self.left_margin+15,self.top_margin+15,"yellow",self.font_size-2,"start","%d CDS" % cds_n))
        return svg
    
    def _draw_infinity_column(self,infinity_ls,shift,reg,dent):
        svg = []
        sign =[-1,1]
        # Add Inf dots
        for i in range(len(infinity_ls)):
            gene,moltype,p,v = infinity_ls[i]
            description = self._format_description(gene,moltype)
            x = 2*self.left_margin+self.plot_width+shift*self.right_margin+random.randint(-10*self.basic_dot_size,10*self.basic_dot_size)
            y = self.top_margin+(.5+self.maxZ-v)*self.side_scale
            r = self.basic_dot_size
            color = self.basic_color
            dm = 0
            if p <= self.p_cutoff:
                dm = 1
            if v >= 2:
                r += 2
                color = self._get_color(moltype,reg,dm)
                self.info[self.categories[reg+2]].append([gene,moltype,p,v*sign[reg]])
                if v >= 4:
                    r += 2
            if gene in self.highlight_genes:
                color = self.highlighted_color 
            strock_color = self.strock_color
            strock_width = self.strock_width
            if gene in self.outlined_genes:
                strock_color = self.highlighted_strock_color
                strock_width += 1
            svg.append(("<a xlink:title=\"%s\">" % description)+"<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"%s\" stroke-width=\"%d\" /></a>" % 
                        (x,y,r,color,strock_color,strock_width))
        if self.categories[reg+2]:
            cds_n = len(list(filter(lambda ls: ls[1].lower() in ("protein_coding","cds") and ls[3] >= 2, self.info[self.categories[reg+2]])))
            rna_n = len(list(filter(lambda ls: ls[1].lower() in ("ncrna","misc_rna","rna") and ls[3] >= 2, self.info[self.categories[reg+2]])))
            if rna_n:
                legend = "%dc+%dnc" % (cds_n,rna_n)
            else:
                legend = "%d CDS" % cds_n
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                (5+2*self.left_margin+self.plot_width+(shift-dent/10.0)*self.right_margin/2+dent-dent,self.top_margin+self.plot_height+15,"yellow",self.font_size-2,"start",legend))
        return svg
    
    def _draw_legend(self):
        svg = []
        # Legend
        svg += self._draw_colorscheme(2*self.left_margin+self.plot_width+self.right_margin/2,self.top_margin+self.plot_height+70)
        return svg

    def _svg_template(self,axisX_title,axisY_title,flg_infinity_plot=True,flg_diagonal=True):
        horizontal_Y = self.height-2*self.bottom_margin
        svg = ["<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" viewbox=\"0 0 %d %d\">" % (self.width,self.height)]
        # background
        if flg_infinity_plot:
            svg.append("<rect x=\"0\" y=\"0\" fill=\"%s\" width=\"%d\" height=\"%d\"/>" % (self.background,self.width,self.height))
        else:
            svg.append("<rect x=\"0\" y=\"0\" fill=\"%s\" width=\"%d\" height=\"%d\"/>" % (self.background,self.width-2*self.right_margin+20,self.height))
        # left main vertical axis
        svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
            (2*self.left_margin,self.top_margin,2*self.left_margin,horizontal_Y+50,"white",2))
        # right main vertical axis
        svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
            (self.width-2*self.right_margin,self.top_margin,self.width-2*self.right_margin,horizontal_Y+50,"white",2))
        if flg_diagonal:
            # slant line
            svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
                (2*self.left_margin,horizontal_Y,self.width-2*self.right_margin,self.top_margin,"white",1))
        if flg_infinity_plot:
            # right RPKM axis
            svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
                (self.width-self.right_margin,self.top_margin,self.width-self.right_margin,horizontal_Y+50,"white",1))
            # horizontal axis
            svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
                (self.left_margin,horizontal_Y,self.width-10,horizontal_Y,"white",2))
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                (self.width-.5*self.right_margin,horizontal_Y+40,"white",self.font_size,"middle",axisY_title))
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                (self.width-1.5*self.right_margin,horizontal_Y+40,"white",self.font_size,"middle",axisX_title))
        else:
            # horizontal axis
            svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
                (self.left_margin,horizontal_Y,self.width-2*self.right_margin,horizontal_Y,"white",2))
        # vertical axis titles
        svg.append(("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" "+
            "transform=\"rotate(%f %d,%d)\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>") %
            (self.left_margin,self.top_margin+(self.height-2*self.bottom_margin-self.top_margin)/2,-90,self.left_margin,self.top_margin+(self.height-2*self.bottom_margin-self.top_margin)/2,
            "white",self.font_size,"middle","%s '%s'" % (self.generic_title,axisY_title)))
        svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
            (self.width-self.right_margin,self.top_margin-1.5*self.font_size,"white",self.font_size,"middle",self.generic_title))
        # horizontal axis titles
        svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
            (self.width-2*self.right_margin-40,horizontal_Y+40,"white",self.font_size,"end","%s '%s'" % (self.generic_title,axisX_title)))

        svg.append("</svg>")
        return svg
    
    def _statistics(self):
        self.x_avr = tools.average(list(map(lambda item: item[-2], self.regulated_genes)))
        self.y_avr = tools.average(list(map(lambda item: item[-1], self.regulated_genes)))
        self._normalize()

        self.maxX = max(list(map(lambda ls: ls[-2], self.regulated_genes)))
        self.minX = min(list(map(lambda ls: ls[-2], self.regulated_genes)))
        self.maxY = max(list(map(lambda ls: ls[-1], self.regulated_genes)))
        self.minY = min(list(map(lambda ls: ls[-1], self.regulated_genes)))
        try:
            self.maxZ = max(list(map(lambda ls: ls[3], self.pos_infinity+self.neg_infinity)))
        except:
            self.maxZ = 0
        try:
            self.minZ = min(list(map(lambda ls: ls[3], self.pos_infinity+self.neg_infinity)))
        except:
            self.minZ = 0
        
        self.x_scale = float(self.plot_width)/(self.maxX-self.minX+1.0)
        self.x_shift = 0.5+abs(self.minX)-abs(int(self.minX))
        self.y_scale = float(self.plot_height)/(self.maxY-self.minY+1.0)
        self.y_shift = 0.5+abs(self.minY)-abs(int(self.minY))
        if self.maxZ or self.minZ:
            self.side_scale = float(self.plot_height)/(self.maxZ-self.minZ+1.0)

        self.pearson = self.expression_degree = 0
        if self.regulated_genes:
            self.pearson = tools.calculate_pearson_correlation(list(map(lambda item: [item[-2],item[-1]], self.regulated_genes)))
        if self.x_avr and self.y_avr:
            self.expression_degree = 45.0-math.atan(self.x_avr/self.y_avr)*180.0/math.pi
        elif self.minX and self.minY:
            self.expression_degree = 45.0 - math.atan(-self.minX/-self.minY)*180.0/math.pi
            
    def _normalize(self):
        return

###############################################################################
class NormExpressionPlot(ExpressionPlot):
    def __init__(self,strain1,strain2,title="",refgene="",filter_file="",highlight_file="",outlined_file="",source_file="",p_cutoff=0.05,infinity_cutoff=0.001):
        self.refgene = refgene
        ExpressionPlot.__init__(self,strain1,strain2,title,filter_file,highlight_file,outlined_file,source_file,p_cutoff,infinity_cutoff)
        self.generic_title = "Log2(baseMean)"

    def set(self,data):
        data = list(filter(lambda item: float(item[2]) and item[6] != "NA", data)) # baseMean != 0
        if self.filter_genes:
            data = list(filter(lambda item: item[0] in self.filter_genes, data))
        self.pos_infinity += list(map(lambda item: [item[0],item[1],float(item[6]),tools.log2(2.0*item[2])], 
            list(filter(lambda item: item[9] < self.infinity_cutoff and item[9] < item[8], data)))) # [gene,moltype,p,2*baseMean]
        self.neg_infinity += list(map(lambda item: [item[0],item[1],float(item[6]),tools.log2(2.0*item[2])], 
            list(filter(lambda item: item[8] < self.infinity_cutoff and item[8] < item[9], data)))) # [gene,moltype,p,2*baseMean]
        used_genes = list(map(lambda item: item[0], self.neg_infinity+self.pos_infinity))
        self.regulated_genes += list(map(lambda item: [item[0],item[1],2**abs(float(item[3])),float(item[6]),
            tools.log2(float(item[8])),tools.log2(float(item[9]))], # [gene,moltype,foldchange,p,RPKM1,RPKM2]
            list(filter(lambda item: item[0] not in used_genes, data))))
            
    def is_reference(self,refgene):
        return list(filter(lambda item: item[0]==refgene, self.regulated_genes))

    def _draw_infinity_column(self,infinity_ls,shift,reg,dent):
        svg = []
        # Add -Inf dots
        for i in range(len(infinity_ls)):
            gene,moltype,v,p = infinity_ls[i]
            description = self._format_description(gene,moltype)
            x = 2*self.left_margin+self.plot_width+shift*self.right_margin+random.randint(-10*self.basic_dot_size,10*self.basic_dot_size)
            y = self.top_margin+(.5+self.maxZ-v)*self.side_scale
            r = self.basic_dot_size
            color = self.basic_color
            dm = 0
            if p <= self.p_cutoff:
                dm = 1
            if abs(v) >= 1:
                r += 2
                color = self._get_color(moltype,reg,dm)
                self.info[self.categories[reg+2]].append(infinity_ls[i])
                if abs(v) >= 2:
                    r += 2
            if gene in self.highlight_genes:
                color = self.highlighted_color 
            strock_color = self.strock_color
            strock_width = self.strock_width
            if gene in self.outlined_genes:
                strock_color = self.highlighted_strock_color
                strock_width += 1
            svg.append(("<a xlink:title=\"%s\">" % description)+"<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"%s\" stroke-width=\"%d\" /></a>" % 
                        (x,y,r,color,strock_color,strock_width))
        if self.categories[reg+2]:
            cds_n = len(list(filter(lambda ls: ls[1].lower() in ("protein_coding","cds"), self.info[self.categories[reg+2]])))
            rna_n = len(list(filter(lambda ls: ls[1].lower() in ("ncrna","misc_rna","rna"), self.info[self.categories[reg+2]])))
            if rna_n:
                legend = "%d+%d" % (cds_n,rna_n)
            else:
                legend = "%d CDS" % cds_n
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                (5+2*self.left_margin+self.plot_width+(shift-dent/10.0)*self.right_margin/2+dent-dent,self.top_margin+self.plot_height+15,"yellow",self.font_size-2,"start",legend))
        return svg
    
    def _draw_trend(self):
        svg = []
        # trend line
        if self.expression_degree:
            if self.minX:
                slop = -self.minX/-self.minY
            else:
                return svg
            if self.expression_degree < 0:
                x = 2*self.left_margin+(self.maxX-self.minX+1.0)*self.x_scale
                y = self.top_margin+(self.maxY-self.maxX/slop)*self.y_scale
            else:
                x = 2*self.left_margin+(self.maxY*slop-self.minX+1.0)*self.x_scale
                y = self.top_margin
            svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
                (2*self.left_margin,self.top_margin+self.plot_height,x,y,"red",1))
        return svg

    def _normalize(self):
        if not self.refgene or not len(list(filter(lambda item: item[0]==self.refgene, self.regulated_genes))):
            self.regulated_genes = list(map(lambda item: item[:4]+[item[4]-self.x_avr,item[5]-self.y_avr], self.regulated_genes))
            self.pos_infinity = list(map(lambda item: item[:3]+[item[3]-self.x_avr], self.pos_infinity))
            self.neg_infinity = list(map(lambda item: item[:3]+[item[3]-self.y_avr], self.neg_infinity))
        else:
            refX,refY = list(filter(lambda item: item[0]==self.refgene, self.regulated_genes))[0][-2:]
            self.regulated_genes = list(map(lambda item: item[:4]+[item[4]-refX,item[5]-refY], self.regulated_genes))
            self.pos_infinity = list(map(lambda item: item[:3]+[item[3]-refX], self.pos_infinity))
            self.neg_infinity = list(map(lambda item: item[:3]+[item[3]-refY], self.neg_infinity))
        self.x_avr = self.y_avr = 0
        
###############################################################################
class BaseMeanPlot(ExpressionPlot):
    def __init__(self,experiment,title="",filter_file="",highlight_file="",outlined_file="",source_file="",p_cutoff=0.05,infinity_cutoff=0.001):
        ExpressionPlot.__init__(self,"","",title,filter_file,highlight_file,outlined_file,source_file,p_cutoff,infinity_cutoff)
        self.experiment = experiment
        self.generic_title = ""

    def set(self,data):
        data = list(filter(lambda item: float(item[2]) and item[3] != "NA" and item[6] != "NA", data)) # baseMean != 0
        if self.filter_genes:
            data = list(filter(lambda item: item[0] in self.filter_genes, data))
        self.regulated_genes += list(map(lambda item: [item[0],item[1],float(item[6]),tools.log2(float(item[2])+1),float(item[3])], # [gene,moltype,p,log2(baseMean),log2(foldchange)]
            data))

    def svg(self): 
        svg = self._svg_template("log2(baseMean+1)","log2(foldChange)",False,False)
        if not self.regulated_genes+self.pos_infinity+self.neg_infinity:
            return svg,""
        
        svg = svg[:-1]
        self._statistics()
        # Draw line scafold and titles
        svg += self._draw_graphs_and_titles()
        # Draw dots
        svg += self._draw_dots()
        # Draw trend line
        svg += self._draw_trend()
        # Draw occurence statistics
        svg += self._draw_occurence_stat()
        # Draw Legends
        svg += self._draw_legend()

        if self.experiment:
            svg.insert(1,"<g id=\"%s\">" % self.experiment)
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s<tspan font-size=\"%f\">%s</tspan></text>\"" % 
                (3*self.left_margin,self.top_margin/2,"white",self.font_size+2,"start",self.experiment,self.font_size,
                " expr. degree = %s" % (tools.format_number(self.expression_degree,1))))
            svg.append("</g>")
        svg.append("</svg>")
        return svg,""
        
        return svg,"\n".join(output)

    def _draw_graphs_and_titles(self):
        svg = []
        if self.title:
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                (self.left_margin+self.plot_width/2,self.height-5,"orange",self.font_size+4,"middle",self.title))
        # X-scale
        for i in range(int(self.maxX-self.minX)+1):
            x = 2*self.left_margin+(self.x_shift+i)*self.x_scale
            if x >= 2*self.left_margin+self.plot_width:
                break
            svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
                (x,self.top_margin+self.plot_height,x,self.top_margin+self.plot_height+5,"white",1))
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%d</text>\"" % 
                (x,self.top_margin+self.plot_height+7+self.font_size,"white",self.font_size,"middle",int(self.minX)+i))
        
        # Y-scale
        for i in range(int(self.maxY-self.minY)+1):
            y = self.top_margin+(.5+self.maxY-int(self.maxY)+i)*self.y_scale
            if y >= self.top_margin + self.plot_height:
                break
            svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
                (2*self.left_margin-5,y,2*self.left_margin,y,"white",1))
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%d</text>\"" % 
                (2*self.left_margin-10,y,"white",self.font_size,"end",int(self.maxY)-i))
        # horizontal line 0
        y = self.top_margin+(.5+self.maxY-self.y_avr)*self.y_scale
        svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
            (2*self.left_margin,y,2*self.left_margin+self.plot_width,y,"grey",1))
        return svg

    def _draw_dots(self):
        svg = []
        # Add dots
        for i in range(len(self.regulated_genes)):
            gene,moltype,p,baseMean,foldchange = self.regulated_genes[i]
            description = self._format_description(gene,moltype)
            x = 2*self.left_margin+(baseMean-self.minX+0.5)*self.x_scale
            y = self.top_margin+(.5+self.maxY-foldchange)*self.y_scale
            r = self.basic_dot_size
            color = self.basic_color
            if abs(foldchange) >= 1:
                reg = dm = 0
                if foldchange > 1:
                    self.info[self.categories[1]].append(self.regulated_genes[i])
                    reg = 1
                elif foldchange < -1:
                    self.info[self.categories[0]].append(self.regulated_genes[i])
                if p <= self.p_cutoff:
                    dm = 1
                color = self._get_color(moltype,reg,dm)
                if abs(foldchange) >= 2:
                    r += 2
                    if abs(foldchange) >= 3:
                        r += 2
            if gene in self.highlight_genes:
                color = self.highlighted_color 
            strock_color = self.strock_color
            strock_width = self.strock_width
            if gene in self.outlined_genes:
                strock_color = self.highlighted_strock_color
                strock_width += 1
            svg.append(("<a xlink:title=\"%s\">" % description)+("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"%s\" stroke-width=\"%d\" /></a>" % 
                        (x,y,r,color,strock_color,strock_width)))
        return svg

    def _draw_trend(self):
        svg = []
        if not len(self.info[self.categories[0]]) or not len(self.info[self.categories[1]]):
            return svg
        # trend line
        X1 = float(sum(list(map(lambda item: item[-2],self.info[self.categories[1]]))))/len(self.info[self.categories[1]])
        Y1 = float(sum(list(map(lambda item: item[-1],self.info[self.categories[1]]))))/len(self.info[self.categories[1]])
        X2 = float(sum(list(map(lambda item: item[-2],self.info[self.categories[0]]))))/len(self.info[self.categories[0]])
        Y2 = float(sum(list(map(lambda item: item[-1],self.info[self.categories[0]]))))/len(self.info[self.categories[0]])
        
        self.expression_degree = 0
        if X1-X2:
            A = (Y1-Y2)/(X1-X2)
            C = ((Y1+Y2) - A*(X1+X2))/2.0
            degree = 180.0*math.atan(abs(A))/math.pi
            if A >= 0:
                self.expression_degree = 90.0-degree
            else:
                self.expression_degree = degree-90.0
        x1 = 2*self.left_margin+(X1-self.minX+0.5)*self.x_scale
        y1 = self.top_margin+(.5+self.maxY-Y1)*self.y_scale
        x2 = 2*self.left_margin+(X2-self.minX+0.5)*self.x_scale
        y2 = self.top_margin+(.5+self.maxY-Y2)*self.y_scale
        x3 = 2*self.left_margin+(-C/A-self.minX+0.5)*self.x_scale

        svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
            (x1,y1,x2,y2,"red",2))
        svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
            (x3,self.top_margin,x3,self.top_margin+self.plot_height,"grey",1))
        return svg
    
    def _draw_legend(self):
        svg = []
        # Legend
        svg += self._draw_colorscheme(2*self.left_margin+self.plot_width-100,self.top_margin+self.plot_height+70)
        return svg

    def _statistics(self):
        self.x_avr = tools.average(list(map(lambda item: item[-2], self.regulated_genes)))
        self.y_avr = tools.average(list(map(lambda item: item[-1], self.regulated_genes)))
        self._normalize()

        self.maxX = max(list(map(lambda ls: ls[-2], self.regulated_genes)))
        self.minX = min(list(map(lambda ls: ls[-2], self.regulated_genes)))
        self.maxY = max(list(map(lambda ls: ls[-1], self.regulated_genes)))
        self.minY = min(list(map(lambda ls: ls[-1], self.regulated_genes)))
        
        self.x_scale = float(self.plot_width)/(self.maxX-self.minX+1.0)
        self.x_shift = 0.5+abs(self.minX)-abs(int(self.minX))
        self.y_scale = float(self.plot_height)/(self.maxY-self.minY+1.0)
        self.y_shift = 0.5+abs(self.minY)-abs(int(self.minY))
        
    def _normalize(self):
        return
            
###############################################################################
if __name__ == "__main__":
    pass
