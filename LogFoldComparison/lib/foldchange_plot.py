import string, random
import tools

global p_val_cutoff
p_val_cutoff = 0.05
global fold_cutoff
fold_cutoff = 1.0
global flg_stroke
flg_stroke = False

def get_stroke_color(p_values):
    if not flg_stroke:
        return "black"
    v1,v2 = p_values
    if v1 == None and v2 == None:
        return "black"
    if ((v1 == None and v2 <= p_val_cutoff) or (v2 == None and v1 <= p_val_cutoff) or
        (v1 > p_val_cutoff and v2 <= p_val_cutoff) or (v2 > p_val_cutoff and v1 <= p_val_cutoff)):
        return "brown"
    if v1 <= p_val_cutoff and v2 <= p_val_cutoff:
        return "red"
    return "black"

def get_color(x_fold,y_fold,basecolor="white"):
    color1 = "orange"               # positive coregulation
    color2 = "moccasin"             # positive Y with no X regulation
    color3 = "plum"                 # alternative regulation positive at Y
    color4 = "olive"                # positive X with no Y regulation
    color5 = "palegreen"            # negative X with no Y regulation
    color6 = "pink"                 # alternative regulation negative at X
    color7 = "powderblue"           # negative Y with no X regulation
    color8 = "mediumpurple"         # negative coregulation
    if x_fold == None:
        if abs(y_fold) < fold_cutoff:
            return basecolor
        if y_fold >= fold_cutoff:
            return color2
        if y_fold <= -fold_cutoff:
            return color7
    if y_fold == None:
        if abs(x_fold) < fold_cutoff:
            return basecolor
        if x_fold >= fold_cutoff:
            return color4
        if x_fold <= -fold_cutoff:
            return color5
    if abs(x_fold) < fold_cutoff and abs(y_fold) < fold_cutoff:
        return basecolor
    if x_fold >= fold_cutoff:
        if y_fold >= fold_cutoff:
            return color1
        if abs(y_fold) < fold_cutoff:
            return color4
        if y_fold <= -fold_cutoff:
            return color3
    if abs(x_fold) < fold_cutoff:
        if y_fold >= fold_cutoff:
            return color2
        if y_fold <= -fold_cutoff:
            return color7
    if x_fold <= -fold_cutoff:
        if y_fold >= fold_cutoff:
            return color6
        if abs(y_fold) < fold_cutoff:
            return color5
        if y_fold <= -fold_cutoff:
            return color8
    
def svg_template(axisX_title,axisY_title,width=800,height=730,left_margin=35,top_margin=50,bottom_margin=130,right_margin=100,background="black",font_size=14):
    horizontal_Y = height-2*bottom_margin
    svg = ["<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" viewbox=\"0 0 %d %d\">" % (width,height)]
    # background
    svg.append("<rect x=\"0\" y=\"0\" fill=\"%s\" width=\"%d\" height=\"%d\"/>" % (background,width,height))
    # left main vertical axis
    svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
        (2*left_margin,top_margin,2*left_margin,horizontal_Y+50,"white",2))
    # right main vertical axis
    svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
        (width-2*right_margin,top_margin,width-2*right_margin,horizontal_Y+50,"white",2))
    # right RPKM axis
    svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
        (width-right_margin,top_margin,width-right_margin,horizontal_Y+50,"white",1))
    # horizontal axis
    svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
        (left_margin,horizontal_Y,width-10,horizontal_Y,"white",2))
    # vertical axis titles
    svg.append(("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" "+
        "transform=\"rotate(%f %d,%d)\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>") %
        (left_margin,top_margin+(height-2*bottom_margin-top_margin)/2,-90,left_margin,top_margin+(height-2*bottom_margin-top_margin)/2,
        "white",font_size,"middle","Log2FoldChange (%s)" % axisY_title))
    svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
        (width-right_margin,top_margin-1.5*font_size,"white",font_size,"middle","Log2FoldChange"))
    # horizontal axis titles
    svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
        (width-2*right_margin-40,horizontal_Y+40,"white",font_size,"end","Log2FoldChange (%s)" % axisX_title))
    svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
        (width-1.5*right_margin,horizontal_Y+40,"white",font_size,"middle",axisX_title))
    svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
        (width-.5*right_margin,horizontal_Y+40,"white",font_size,"middle",axisY_title))

    svg.append("</svg>")
    return svg,left_margin,top_margin,bottom_margin,right_margin,width-2*right_margin-2*left_margin,height-2*bottom_margin-top_margin,font_size

def plot(values,plot_title="",    # values = [[[fold1,p1,goterm1,title1],[fold2,p2,goterm2,title2]...], [[-inf1_RPKM,goterm1,title1],[-inf2_RPKM,goterm2,title2],...], [[+inf1_RPKM,goterm1,title1],...] 
    p_cutoff="No",basic_color = "white",basic_dot_size = 3):
    
    try:
        strain1,strain2 = list(map(lambda s: s.strip(), plot_title.split("|")))
    except:
        strain1,strain2 = list(map(lambda s: s.strip(), plot_title.split("_vs_")))
    global p_val_cutoff
    if p_cutoff != "No":
        p_val_cutoff = float(p_cutoff)
    else:
        p_val_cutoff = 0.05
    categories = ["Co-regulation","Counter-regulation","Expressed in %s" % strain1,"Expressed in %s" % strain2]
    info = dict(zip(categories,list(map(lambda i: [], range(len(categories))))))
    
    svg,left_margin,top_margin,bottom_margin,right_margin,width,height,font_size = svg_template(strain1,strain2)    
    svg = svg[:-1]
    values[0] = list(map(lambda ls: [float(ls[0]),float(ls[1]),ls[2],ls[3],ls[5],ls[4]], values[0]))
    if values[1]:
        values[1] = list(map(lambda ls: [float(ls[0]),ls[1],ls[2],ls[4],ls[3]], values[1]))
    if values[2]:
        values[2] = list(map(lambda ls: [float(ls[0]),ls[1],ls[2],ls[4],ls[3]], values[2]))
    maxX = minX = maxY = minY = 0
    if values[0]:
        maxX = max(list(map(lambda ls: ls[0], values[0])))
        minX = min(list(map(lambda ls: ls[0], values[0])))
        maxY = max(list(map(lambda ls: ls[1], values[0])))
        minY = min(list(map(lambda ls: ls[1], values[0])))
    try:
        max_foldchange = max(list(map(lambda ls: ls[0], values[1]+values[2])))
    except:
        max_foldchange = 0
    try:
        min_foldchange = min(list(map(lambda ls: ls[0], values[1]+values[2])))
    except:
        min_foldchange = 0
    # X-scale
    x_scale = float(width)/(maxX-minX+1.0)
    x_shift = 0.5+abs(minX)-abs(int(minX))
    for i in range(int(maxX-minX)):
        x = 2*left_margin+(x_shift+i)*x_scale
        svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
            (x,top_margin+height,x,top_margin+height+5,"white",1))
        svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%d</text>\"" % 
            (x,top_margin+height+7+font_size,"white",font_size,"middle",int(minX)+i))
    # vertical line -1
    i = -fold_cutoff-int(minX)
    x = 2*left_margin+(x_shift+i)*x_scale
    x_vline_minus = x
    svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
        (x,top_margin,x,top_margin+height+5,"grey",1))
    # vertical line 1
    i = fold_cutoff-int(minX)
    x = 2*left_margin+(x_shift+i)*x_scale
    x_vline_plus = x
    svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
        (x,top_margin,x,top_margin+height+5,"grey",1))
    
    # Y-scale
    y_scale = float(height)/(maxY-minY+1.0)
    y_shift = 0.5+abs(minY)-abs(int(minY))
    for i in range(int(maxY-minY)+1):
        y = top_margin+(.5+maxY-int(maxY)+i)*y_scale
        svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
            (2*left_margin-5,y,2*left_margin,y,"white",1))
        svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%d</text>\"" % 
            (2*left_margin-10,y,"white",font_size,"end",int(maxY)-i))
    # horizontal line -1
    y = top_margin+(fold_cutoff+.5+maxY)*y_scale
    y_hline_minus = y
    svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
        (2*left_margin,y,2*left_margin+width,y,"grey",1))
    # horizontal line 1
    y = top_margin+(maxY-fold_cutoff+.5)*y_scale
    y_hline_plus = y
    svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
        (2*left_margin,y,2*left_margin+width,y,"grey",1))
        
    # Side-scale
    if max_foldchange+min_foldchange:
        side_scale = float(height)/(max_foldchange-min_foldchange+1.0)
        side_axis_shift = 0.5+abs(min_foldchange)-abs(int(min_foldchange))
        # side horizontal line -1
        y = top_margin+(1.5+max_foldchange)*side_scale
        svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
            (2*left_margin+width,y,2*left_margin+width+2*right_margin,y,"grey",1))
        # side horizontal line 1
        y = top_margin+(max_foldchange-.5)*side_scale
        svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
            (2*left_margin+width,y,2*left_margin+width+2*right_margin,y,"grey",1))
        for i in range(int(max_foldchange-min_foldchange)):
            #y = top_margin+(side_axis_shift+i-.5)*side_scale
            y = top_margin+(.5+max_foldchange-int(max_foldchange)+i)*side_scale
            svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%d\" />" % 
                (2*left_margin+width+right_margin-5,y,2*left_margin+width+right_margin+5,y,"white",1))
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%d</text>\"" % 
                (2*left_margin+width+right_margin+10,y,"white",font_size,"start",int(max_foldchange)-i))
    
    # Add volcano dots
    dots = [[],[],[],[]]
    for i in range(len(values[0])):
        k = 3
        fold1,fold2,gene,tag,p_values,moltype = values[0][i]
        title = "%s; %s" % (tag,gene)
        x = 2*left_margin+(fold1-minX+0.5)*x_scale
        y = top_margin+(.5+maxY-fold2)*y_scale
        r = basic_dot_size
        if abs(fold1) < fold_cutoff and abs(fold2) < fold_cutoff:
            color = basic_color
        else:
            color = get_color(fold1,fold2)
            r += 1
            k -= 1
            if fold1 < 0 and fold2 > 0 or fold1 > 0 and fold2 < 0:
                info[categories[1]].append(values[0][i])
            elif fold1 < 0 and fold2 < 0 or fold1 > 0 and fold2 > 0:
                info[categories[0]].append(values[0][i])
            if p_values[0] != None and p_values[0] <= p_val_cutoff:
                r += 1
                k -= 1
            if p_values[1] != None and p_values[1] <= p_val_cutoff:
                r += 2
                k -= 1
        if color == basic_color:
            stroke_color = "black"
        else:
            stroke_color = get_stroke_color(p_values)
        stroke_width = 1
        if stroke_color != "black":
            stroke_width = 2
        dots[k].append("<a xlink:title=\"%s\"><circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"%s\" stroke-width=\"%d\" /></a>" % 
                    (title,x,y,r,color,stroke_color,stroke_width))
    for k in range(4):
        svg += dots[k]

    # Stat values
    contingency_table = tools.contingency_table(info[categories[0]]+info[categories[1]],fold_cutoff)
    if contingency_table[0]:
        if contingency_table[0][1]:
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                (2*left_margin+width-15,top_margin+5,"yellow",font_size-2,"end","%d+%d" % (contingency_table[0][0],contingency_table[0][1])))
        else:
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                (2*left_margin+width-15,top_margin+5,"yellow",font_size-2,"end","%d CDS" % contingency_table[0][0]))
    if contingency_table[1]:
        if contingency_table[1][1]:
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                (2*left_margin+15,top_margin+5,"yellow",font_size-2,"start","%d+%d" % (contingency_table[1][0],contingency_table[1][1])))
        else:
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                (2*left_margin+15,top_margin+5,"yellow",font_size-2,"start","%d CDS" % (contingency_table[1][0])))
    if contingency_table[2]:
        if contingency_table[2][1]:
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                (2*left_margin+width-15,y_hline_minus-5,"yellow",font_size-2,"end","%d+%d" % (contingency_table[2][0],contingency_table[2][1])))
        else:
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                (2*left_margin+width-15,y_hline_minus-5,"yellow",font_size-2,"end","%d CDS" % (contingency_table[2][0])))
    if contingency_table[3]:
        if contingency_table[3][1]:
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                (2*left_margin+15,y_hline_minus-5,"yellow",font_size-2,"start","%d+%d" % (contingency_table[3][0],contingency_table[3][1])))
        else:
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                (2*left_margin+15,y_hline_minus-5,"yellow",font_size-2,"start","%d CDS" % (contingency_table[3][0])))
    if contingency_table[4]:
        if contingency_table[4][1]:
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                (x_vline_minus+5,top_margin+5,"yellow",font_size-2,"start","%d+%d" % (contingency_table[4][0],contingency_table[4][1])))
        else:
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                (x_vline_minus+5,top_margin+5,"yellow",font_size-2,"start","%d CDS" % (contingency_table[4][0])))
    if contingency_table[5]:
        if contingency_table[5][1]:
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                (x_vline_minus+5,top_margin+height-5,"yellow",font_size-2,"start","%d+%d" % (contingency_table[5][0],contingency_table[5][1])))
        else:
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                (x_vline_minus+5,top_margin+height-5,"yellow",font_size-2,"start","%d CDS" % (contingency_table[5][0])))
    if contingency_table[7]:
        if contingency_table[7][1]:
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                (2*left_margin+width-15,top_margin+height-5,"yellow",font_size-2,"end","%d+%d" % (contingency_table[7][0],contingency_table[7][1])))
        else:
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                (2*left_margin+width-15,top_margin+height-5,"yellow",font_size-2,"end","%d CDS" % (contingency_table[7][0])))
    if contingency_table[8]:
        if contingency_table[8][1]:
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                (2*left_margin+15,top_margin+height-5,"yellow",font_size-2,"start","%d+%d" % (contingency_table[8][0],contingency_table[3][1])))
        else:
            svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
                (2*left_margin+15,top_margin+height-5,"yellow",font_size-2,"start","%d CDS" % contingency_table[8][0]))
    total_cds = sum(list(map(lambda i: contingency_table[i][0], range(len(contingency_table)))))
    total_ncrna = sum(list(map(lambda i: contingency_table[i][1], range(len(contingency_table)))))
    if total_ncrna:
        svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
            (5,top_margin/2,"yellow",font_size,"start","%d+%d" % (total_cds,total_ncrna)))
    else:
        svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
            (5,top_margin/2,"yellow",font_size,"start","%d CDS" % total_cds))

    # Add -Inf dots
    dots = [[],[],[],[]]
    for i in range(len(values[1])):
        k = 3
        v,gene,tag,p_values,moltype = values[1][i]
        title = "%s; %s" % (tag,gene)
        x = 2*left_margin+width+right_margin/2+random.randint(-10*basic_dot_size,10*basic_dot_size)
        y = top_margin+(.5+max_foldchange-v)*side_scale
        r = basic_dot_size
        color = basic_color
        if abs(v) >= fold_cutoff:
            r += 1
            k -= 1
            color = get_color(v,None)
            info[categories[2]].append(values[1][i])
        if p_values[0] != None and p_values[0] <= p_val_cutoff:
            r += 1
            k -= 1
        if p_values[1] != None and p_values[1] <= p_val_cutoff:
            r += 2
            k -= 1
        if color == basic_color:
            stroke_color = "black"
        else:
            stroke_color = get_stroke_color(p_values)
        stroke_width = 1
        if stroke_color != "black":
            stroke_width = 2
        dots[k].append("<a xlink:title=\"%s\"><circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"%s\" stroke-width=\"%d\" /></a>" % 
                    (title,x,y,r,color,stroke_color,stroke_width))
    for k in range(4):
        svg += dots[k]

    if values[1]:
        cds_up = len(list(filter(lambda ls: ls[0] >= 1 and ls[-1].upper()=="CDS", values[1])))
        cds_down = len(list(filter(lambda ls: ls[0] <= -1 and ls[-1].upper()=="CDS", values[1])))
        rna_up = len(list(filter(lambda ls: ls[0] >= 1 and ls[-1].upper() in ("RNA","NCRNA"), values[1])))
        rna_down = len(list(filter(lambda ls: ls[0] <= -1 and ls[-1].upper() in ("RNA","NCRNA"), values[1])))
        if rna_up:
            legend = "%d+%d" % (cds_up,rna_up)
        else:
            legend = "%d CDS" % cds_up
        svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
            (2*left_margin+width+5,top_margin+5,"yellow",font_size-2,"start",legend))
        if rna_down:
            legend = "%d+%d" % (cds_down,rna_down)
        else:
            legend = "%d CDS" % cds_down
        svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
            (2*left_margin+width+5,top_margin+height+15,"yellow",font_size-2,"start",legend))
            
    # Add +Inf dots
    dots = [[],[],[],[]]
    for i in range(len(values[2])):
        k = 3
        v,gene,tag,p_values,moltype = values[2][i]
        title = "%s; %s" % (tag,gene)
        x = 2*left_margin+width+right_margin/2+right_margin+random.randint(-10*basic_dot_size,10*basic_dot_size)
        y = top_margin+(.5+max_foldchange-v)*side_scale
        r = basic_dot_size
        color = basic_color
        if abs(v) >= fold_cutoff:
            r += 1
            k -= 1
            color = get_color(None,v)
            info[categories[3]].append(values[2][i])
        if p_values[0] != None and p_values[0] <= p_val_cutoff:
            r += 1
            k -= 1
        if p_values[1] != None and p_values[1] <= p_val_cutoff:
            r += 2
            k -= 1
        if color == basic_color:
            stroke_color = "black"
        else:
            stroke_color = get_stroke_color(p_values)
        stroke_width = 1
        if stroke_color != "black":
            stroke_width = 2
        svg.append("<a xlink:title=\"%s\"><circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"%s\" stroke-width=\"%d\" /></a>" % 
                    (title,x,y,r,color,stroke_color,stroke_width))
        
    for k in range(4):
        svg += dots[k]

    if values[2]:
        cds_up = len(list(filter(lambda ls: ls[0] >= 1 and ls[-1].upper()=="CDS", values[2])))
        cds_down = len(list(filter(lambda ls: ls[0] <= -1 and ls[-1].upper()=="CDS", values[2])))
        rna_up = len(list(filter(lambda ls: ls[0] >= 1 and ls[-1].upper() in ("RNA","NCRNA"), values[2])))
        rna_down = len(list(filter(lambda ls: ls[0] <= -1 and ls[-1].upper() in ("RNA","NCRNA"), values[2])))
        if rna_up:
            legend = "%d+%d" % (cds_up,rna_up)
        else:
            legend = "%d CDS" % cds_up
        svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
            (2*left_margin+width+2*right_margin-5,top_margin+5,"yellow",font_size-2,"end",legend))
        if rna_down:
            legend = "%d+%d" % (cds_down,rna_down)
        else:
            legend = "%d CDS" % cds_down
        svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" % 
            (2*left_margin+width+2*right_margin-5,top_margin+height+15,"yellow",font_size-2,"end",legend))
    
    # Legend
    x = 3*left_margin
    y = top_margin+height+75
    y_step = 20
    x_shift = 0
    r = basic_dot_size+4
    
    color1 = "orange"               # positive coregulation
    color2 = "moccasin"             # positive Y with no X regulation
    color3 = "plum"                 # alternative regulation positive at Y
    color4 = "olive"                # positive X with no Y regulation
    color5 = "palegreen"            # negative X with no Y regulation
    color6 = "pink"                 # alternative regulation negative at X
    color7 = "powderblue"           # negative Y with no X regulation
    color8 = "mediumpurple"         # negative coregulation
    ls_goterms = [["positive coregulation",fold_cutoff+1,fold_cutoff+1],
        ["upregulation in %s and downregulation in %s" % (strain1,strain2),fold_cutoff+1,0],
        ["upregulation in %s and downregulation in %s" % (strain1,strain2),fold_cutoff+1,-fold_cutoff-1],
        ["upregulation in %s and downregulation in %s" % (strain2,strain1),0,fold_cutoff+1],
        ["upregulation in %s and downregulation in %s" % (strain2,strain1),-fold_cutoff-1,fold_cutoff+1],
        ["downregulation in %s and no regulation in %s" % (strain1,strain2),-fold_cutoff-1,0],
        ["downregulation in %s and no regulation in %s" % (strain2,strain1),0,-fold_cutoff-1],
        ["negative coregulation",-fold_cutoff-1,-fold_cutoff-1]]
    for i in range(len(ls_goterms)):
        if i == 8:
            y = top_margin+height+75
            x_shift = 300
        goterm,a,b = ls_goterms[i]
        svg.append("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"%s\" stroke-width=\"1\" />" % 
                    (x+x_shift,y,r,get_color(a,b),"black"))
        svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>\"" %
            (x+x_shift+15,y,"white",font_size-2,"start","- %s" % goterm))
        y += y_step
    
    print(tools.calculate_pearson_correlation(list(map(lambda item: [item[0],item[1]], info[categories[0]]+info[categories[1]]))))
    #distance = tools.calculate_distance(list(map(lambda item: [item[0],item[1]], info[categories[0]]+info[categories[1]])))
    pearson,p_value = tools.calculate_pearson_correlation(list(map(lambda item: [item[0],item[1]], info[categories[0]]+info[categories[1]])))
    
    if plot_title:
        svg.insert(1,"<g id=\"%s\">" % plot_title)
        svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"%s\" font-size=\"%d\" style=\"text-anchor:%s\">%s<tspan font-size=\"%f\">%s</tspan></text>\"" % 
            (4*left_margin,top_margin/2,"white",font_size+2,"start",plot_title,font_size,
            " Pearson Corr. = %s (p-value = %s)" % (tools.format_number(pearson,3),tools.format_number(p_value,4))))
        svg.append("</g>")
    svg.append("</svg>")
    output = []
    for i in range(2):
        key = categories[i]
        info[key].sort(key=lambda ls: ls[0]**2+ls[1]**2, reverse=True)
        output.append("%s\n\t%s" % (key,"\n\t".join(list(map(lambda item: "\t".join(list(map(lambda v: str(v), item))), info[key])))))
    for i in range(2,4,1):
        key = categories[i]
        info[key].sort(key=lambda ls: ls[0], reverse=True)
        output.append("%s\n\t\t%s" % (key,"\n\t\t".join(list(map(lambda item: "\t".join(list(map(lambda v: str(v), item))), info[key])))))
    
    return "\n".join(svg),"\n".join(output)

###############################################################################
if __name__ == "__main__":
    import seq_io
    IO = seq_io.IO()
    #svg,left_margin,top_margin,right_margin,width,height,font_size = svg_template()
    svg = volcano([[[-7.33,0.0003],[6.8,0.000001],[-2.33,0.32],[4.4,0.04],[3.7,0.033]],[123,178,2005,18],[91,118,576,23]],"Test")
    IO.save("\n".join(svg),"volcano_plot.svg")
