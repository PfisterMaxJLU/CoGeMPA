import argparse
import json
import re
import matplotlib.pyplot as plt
import numpy as np
import csv
import re
import requests
import xlsxwriter
from appJar import gui
import subprocess
from datetime import datetime
import os

#use Objects to carry data next time

def argp():
    parser = argparse.ArgumentParser(description='Input for various functions of this program')

    group= parser.add_mutually_exclusive_group()
    group.add_argument('-a', '--met_analysis', action='store_true', help= "For analysing the KofamKOALA output (-dbj [first time and when it should be changed] -f1 required)")
    group.add_argument('-c', '--csv_analysis', action='store_true', help= "For analysing an EDGAR csv file (-csv required)")
    group.add_argument('-r', '--rem_duplicates', action='store_true', help= "For removing duplicates in the EDGAR faa. file (-e -o required) [or any other faa. file]")
    group.add_argument('-d', '--make_database', action='store_true', help= "For creating the json \"database\" for this program (-db required)")
    group.add_argument('-i', '--install', action='store_true', help= "Auto install KEGG DB")  
    group.add_argument('-cm', '--csv_meta', action='store_true', help= "Metabolic analysis of csv file (-dbj [first time and when it should be changed] -csv and -f1 required)")
    group.add_argument('-rko', '--remko', action='store_true', help= "Removes duplicate KO numbers in KofamScan output -f1 required, output name set automatic. Important Info: Also removes Gen \"duplicates\" within a single organism's set, if K number is deemed to be identical! ")  

    group.add_argument('-sg', '--start_gui', action='store_true', help= "For starting the GUI version of this program")

    parser.add_argument('-exl', '--export_excel', action='store_true', help= "Exports a sorted list of gen descriptions as an xls file for -a and -cm (-o can! be used)") 

    parser.add_argument('-fa', '--fullauto', action='store_true', help= "The selected files will be first processed in KOFAM-Scan then analysed (-f1 required")
    parser.add_argument('-ra', '--rem_auto', action='store_true', help= "Use if there are duplicate EDGAR-Id's jamming auto mode") 

    parser.add_argument('-e', '--EDGAR', metavar='', help= "(EDGAR) faa. file where duplicates will be removed")
    parser.add_argument('-o', '--output', metavar='', help= "Name that the output file will have")
    parser.add_argument('-db', '--database', metavar='', help= "Path to original KEGG \"database\" (Download: https://www.kegg.jp/kegg-bin/get_htext?ko00001.keg), or the change brite type for -i")
    
    parser.add_argument('-dbj', '--dbjson', metavar='', default=None, help= "Name for converted KEGG \"database\" json.")

    parser.add_argument('-csv', '--csvpath', metavar='', help= "Path to EDGAR csv which will be analysed")

    parser.add_argument('-f1', '--file1', metavar='', help= "Path to KofamKOALA which will be analysed")
    parser.add_argument('-f2', '--file2', metavar='', help= "Path to Optional second file")
    parser.add_argument('-f3', '--file3', metavar='', help= "Path to Optional third file")

    parser.add_argument('-na','--not_annotated', action='store_true', help= "Set this if you want to see the fraction of KEGG sequences not annotated in the KEGG database")
    parser.add_argument('-pie','--pie_chart', action='store_true', help= "Set this if you want pie charts instead of bar charts")
    parser.add_argument('-per', '--percent', action='store_true', help= "Set this if you want Bar Charts to show relative percent in their Group/Species instead of total numbers")
    parser.add_argument('-jf', '--just_first', action='store_true', help= "Use this if you want the program to just the first \"most relevant\" categories for a single KA-number")  

    #f1 for all maybe?

    args = parser.parse_args()
    return(args)


def analyse(KEGG_db_path, EDGAR_paths, set_not_annotated, pie, output_path, wants_xlsx, percent, just_first, gui_on=False, app=None): #EDGAR_path also known as KOFAM_file_path not consistent, also "gui_on" now redundant because of app

    detail_list = {}

    try:
        EDGAR_paths = list(filter(None, EDGAR_paths)) 

        plt.close('all')

        with open(KEGG_db_path, 'r') as db:
            database = json.load(db)

        with open("./categories_configuration.json", 'r') as config: 
            cat_cfg = json.load(config)


        for i, p in enumerate(EDGAR_paths):

            with open(p, "r") as EDGAR_file:
        
                c = count_categorie(EDGAR_file, database, cat_cfg["names"], p, detail_list, just_first)

                plotter(c, cat_cfg["names"], cat_cfg["colors"], i, len(EDGAR_paths), p, set_not_annotated, pie, gui_on, detail_list, None, app, percent)
            
        if pie == False:
            legend_properties = {'weight':'bold', 'size':14}
            plt.legend(prop=legend_properties)

        # with open('sorted_gen_descriptions', 'w') as fp: #here for debugging
        #     json.dump(detail_list, fp, indent=4)

        if wants_xlsx:
            xlsx_export(detail_list, cat_cfg["names"], output_path)

        plt.show()


    except Exception as e: 
        print(e)
        print("Please select the appropriate Database or files") 


def plotter(data_dict, order, colors, num, maxn, filename, set_not_annotated, pie, gui_on, detail_list, KOFAM_file_path, app, percent):
    '''
    foo bar
    -------
    '''
    
    values = []

    for key in order:
        values.append(data_dict[key])
        
    if set_not_annotated == False:
        values = values[:-1]
    if set_not_annotated == False:
        order = order[:-1]
    if set_not_annotated == False:
        colors = colors[:-1]
    
    #bad use name
    
    if pie == False:
        
        if percent:
            values_total = values 
            values = values.copy()
            total_to_percent(values)
            plt.ylabel('Relative (%)', fontweight='bold', fontsize=18)
        else:
            plt.ylabel('Amount', fontweight='bold', fontsize=18)

        x_pos = np.arange(len(values))*maxn*1
        width = 0.8
        fig = plt.gcf()
        p = plt.bar(x_pos + width*num , values, align='center', alpha=1-num*0.65/(maxn), color=colors, width=width)
        
        cut_EDGAR_path = os.path.basename(filename)

        if not percent:
            p.set_label(f"{cut_EDGAR_path} ---> Total found categories: {sum(values)}")
            
        else:
            p.set_label(f"{cut_EDGAR_path} ---> Total found categories: {sum(values_total)}")


        annot = plt.annotate(text="", xy=(0,0), xytext=(-20,20),textcoords="offset points",
                    bbox=dict(boxstyle="round", fc="black", ec="b", lw=2),
                    arrowprops=dict(arrowstyle="->"))
        annot.set_visible(False)

        plt.yticks(weight='bold', rotation=90, verticalalignment='center', fontsize=18)

        #plt.tight_layout()
        plt.gcf().subplots_adjust(top=0.975)
        plt.gcf().subplots_adjust(bottom=0.53)
        plt.gcf().subplots_adjust(left=0.03)
        plt.gcf().subplots_adjust(right=0.99)
        
        plt.xticks(x_pos + width*num*0.5, order, rotation=90, fontweight='bold', fontsize=18)
        

        def update_annot(bar, i):
            x = bar.get_x()+bar.get_width()/2.
            y = bar.get_y()+bar.get_height()
            annot.xy = (x,y)
            if percent:
                text = f"({y}% ≙ {values_total[i]})"
            else:
                text = f"({y})"
            annot.set_text(text)
            annot.get_bbox_patch().set_alpha(0.4)


        def on_hover(event):
            vis = annot.get_visible()
            for i, bar in enumerate(p):
                cont, _ = bar.contains(event)
                if cont:
                    update_annot(bar, i)
                    annot.set_visible(True)
                    fig.canvas.draw_idle()
                    return
            if vis:
                annot.set_visible(False)
                fig.canvas.draw_idle()


        def on_click(event):
            close = False
            for i, bar in enumerate(p):
                cont, _ = bar.contains(event)
                if cont:
                    if app != None:
                        app.setLabel("infolabel", f"{order[i]}, {filename}, ({KOFAM_file_path}) | data exported")
                        if order[i] == 'Not Annotated':
                            app.setLabel("infolabel", "It is not possible to export Gen-descriptions for Genes that were not annotated")
                    onclick_xlsx_export(order[i], num, filename, gui_on, detail_list, KOFAM_file_path, close)

        #add legend picker maybe
        # def on_close(event): rip maybe in v2
        #     try:
        #         onclick_xlsx_export(workbook, None, None, None, None, None, None, True)
        #         print("xlsx saved")
        #     except Exception as e: 
        #         print(e)

        fig.canvas.mpl_connect("motion_notify_event", on_hover)
        fig.canvas.mpl_connect("button_press_event", on_click)
        #fig.canvas.mpl_connect('close_event', on_close)


    if pie == True: 
        values, order, colors = value_order_colour_rearranger(values, order, colors)
        #ohboy

        values.reverse() 
        order.reverse()
        colors.reverse()

        counter = 0

        if sum(values) == 0:
            values = [1]
            order = ["No Genes in this Category"]
            colors = ["#000000"]
        
        for n in values:
            if n > 0:
                counter += 1

        fig, (ax1) = plt.subplots(sharey=True)

        for i, o in enumerate(values):
            order[i] = f"{order[i]} --> ({o})"
        
        ax1.pie(values, autopct=percent_maker, startangle=90, colors=colors, pctdistance=0.7, counterclock=False, wedgeprops={"edgecolor":"k",'linewidth': 0.35, 'antialiased': True})

        plt.gcf().subplots_adjust(right=0.56)

        legend_properties = {'weight':'bold', 'size':14}
        ax1.legend(order[:counter], bbox_to_anchor=(0.95, 0.8), loc='upper left', borderaxespad=0.1, prop=legend_properties)

        if KOFAM_file_path == None:
            nameospath = os.path.basename(filename)
            plt.title(f'Metabolic Profile of:\n{nameospath}\n Total found categories: {sum(values)}', fontweight='bold', fontsize=14)
            
        else:
            nameospath = os.path.basename(KOFAM_file_path)
            plt.title(f'Metabolic Profile from the {filename} of:\n{nameospath}\n Total found categories: {sum(values)}', fontweight='bold', fontsize=14)      

def total_to_percent(values): #relativ to their group! species! #rip pefekt 40kb because of this fix
    val_sum = sum(values)
    if val_sum != 0:
        for i, v in enumerate(values):
            new_val = (v/val_sum)*100
            new_val = round(new_val, 2)
            values[i] = new_val


def percent_maker(val):

    a  = np.round(val, 1)

    if a >= 2:
        percent = f"{a}%"
    else:
        percent = ""

    return percent
    

def value_order_colour_rearranger(values, order, colors): 
    
    combined = zip(values, order, colors)
    combined = sorted(combined, key=lambda x: x[0])

    values = [x[0] for x in combined]
    order = [x[1] for x in combined]
    colors = [x[2] for x in combined]


    # sorted_values = sorted(values) #jank recursive bubble sort

    # for i, e in enumerate(values):
    #     if e < values[i-1]:
    #         if i-1 >= 0:
    #             values[i], values[i-1] =  values[i-1], values[i]
    #             order[i], order[i-1] =  order[i-1], order[i]
    #             colors[i], colors[i-1] =  colors[i-1], colors[i]

    # if values != sorted_values:
    #     value_order_colour_rearranger(values, order, colors)

    return values, order, colors
    

def count_categorie(EDGAR_file, database, categories, KOFAM_file_path, detail_list, just_first):

    regexp = re.compile(r'[\t][K]{1}[0-9]{5}')
    not_annotated = 0
    not_in_db = 0
    found_K_categories = {}

    for c in categories:
        found_K_categories[c] = 0

    #make found_K_categories with database categories, sometimes not all there

    for line in EDGAR_file:
        #get edgar ids here

        splitline = line.split("\t")
                    
        EDGAR_ID = splitline[0].strip("\n")


        if regexp.search(line):
            current_K_number = line[-7:-1]

            if current_K_number in database:

                if just_first:

                    if database[current_K_number]['categorie_sub'][0] in categories:
                        found_K_categories[database[current_K_number]['categorie_sub'][0]] += 1
                        save_details(KOFAM_file_path, 3, database[current_K_number]['categorie_sub'][0], database, current_K_number, detail_list, EDGAR_ID)

                    else: 
                            found_K_categories["Other"] += 1
                            save_details(KOFAM_file_path, 3, 'Other', database, current_K_number, detail_list, EDGAR_ID)

                else:

                    for c in database[current_K_number]['categorie_sub']:
                        
                        
                        if c in categories:
                            found_K_categories[c] += 1
                            save_details(KOFAM_file_path, 3, c, database, current_K_number, detail_list, EDGAR_ID)

                        else: 
                            found_K_categories["Other"] += 1
                            save_details(KOFAM_file_path, 3, 'Other', database, current_K_number, detail_list, EDGAR_ID)

            else:
                not_in_db += 1 # tecnicaly not false, i mean if there is a K number it is something, here could bee a in seconday db search
                #print("ERROR K-number not in database") 
        else:
            not_annotated += 1 #also return this at some point, show in bar chart
    
    
    found_K_categories["Not Annotated"] = not_annotated

    print(f"{not_in_db} sequences were not found in your local DB")
    print(f"{not_annotated} sequences were not annotated in the KofamKOALA output")
    
    return found_K_categories


def pan_core_analysis(csv_path, pie):

    try:
        core_gen = 0
        disp_gen = 0
        sing_gen = 0


        with open(csv_path, "r") as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            for lines in csv_reader:
                len_csv = len(lines) - 2

                missing_gen = 0
                missing_gen = lines.count(' - , -')

                if len_csv - missing_gen == len_csv:
                    core_gen += 1
                elif len_csv - missing_gen == 1:
                    sing_gen += 1
                else:
                    disp_gen += 1

        values_plot = [core_gen,disp_gen,sing_gen]

        x_pos = np.arange(len(values_plot))

        plt.close('all')

        colors = ["#99bb55","#ee9944","#444466"]
        description = ["Core Gens", "Dispensable Gens", "Singleton Gens"]
        
        if pie == False:
            p = plt.bar(x_pos, values_plot, align='center', color=colors)
            plt.xticks(x_pos, description, fontweight='bold', fontsize=18)
            plt.yticks(fontweight='bold', fontsize=18)
            plt.ylabel('Amount', fontweight='bold', fontsize=18)
            # no time to add interactivity here, maybe later, planing!
        
        elif pie == True:

            for i, o in enumerate(values_plot):
                description[i] = f"{description[i]} --> ({o})"

            fig, (ax1) = plt.subplots(sharey=True)
            ax1.pie(values_plot , autopct=percent_maker, startangle=90, colors=colors, pctdistance=0.7, counterclock=False, wedgeprops={"edgecolor":"k",'linewidth': 0.35, 'antialiased': True})
            legend_properties = {'weight':'bold', 'size':14}
            ax1.legend(description, bbox_to_anchor=(0.95, 0.75), loc='upper left', borderaxespad=0.1, prop=legend_properties)
            nameospath = os.path.basename(csv_path)
            plt.title(f'Distribution of Core-,Disp.- and Sing.- Genes in:\n{nameospath} \n Total found categories: {sum(values_plot)}', fontweight='bold', fontsize=14)

        plt.show()
    except Exception as e: 
        print(e)
        print("Please select an appropriate .csv")


def pan_core_analysis_meta(csv_path, KEGG_db_path, KOFAM_file_path, set_not_annotated, pie, output_path, wants_xlsx, percent, just_first, gui_on=False, app=None): # try catch!

    try:
        singelton_dic = {}
        dispensable_dic = {}
        core_dic = {}
        detail_list = {}

        singelton_dic = dic_initialiser(singelton_dic)
        dispensable_dic = dic_initialiser(dispensable_dic)
        core_dic = dic_initialiser(core_dic)


        with open(csv_path, "r") as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            
            group_dic = group_lister(csv_reader)

        with open("./categories_configuration.json", 'r') as config: 
            cat_cfg = json.load(config)

            with open(KEGG_db_path, 'r') as db:
                database = json.load(db)

                with open(KOFAM_file_path, "r") as KOFAM_file: # merge them later

                    
                    for line in KOFAM_file:

                        splitline = line.split("\t")
                        
                        if len(splitline) == 1:
                            EDGAR_ID = splitline[0].strip("\n")

                            dict_merger(EDGAR_ID, None, group_dic, database, singelton_dic, dispensable_dic, core_dic, KOFAM_file_path, detail_list, just_first)

                        else:
                            EDGAR_ID = splitline[0].strip("\n")
                            KEGG_ID = splitline[1].strip("\n")
                            
                            dict_merger(EDGAR_ID, KEGG_ID, group_dic, database, singelton_dic, dispensable_dic, core_dic, KOFAM_file_path, detail_list, just_first)

            dic_list = [core_dic, dispensable_dic, singelton_dic]

            plt.close('all')
            taglist = iter(["Core Gens", "Dispensable Gens", "Singelton Gens"]) #Ja hier dreht es sich...für den plot
            
            
            for i, dic in enumerate(dic_list):
                tag =  next(taglist)

                plotter(dic, cat_cfg["names"], cat_cfg["colors"], i, len(dic_list), tag, set_not_annotated, pie, gui_on, detail_list, KOFAM_file_path, app, percent) #here set name

        if pie == False:
            legend_properties = {'weight':'bold', 'size':14}
            plt.legend(prop=legend_properties)
    except Exception as e: 
        print(e)

    # with open('sorted_gen_descriptions', 'w') as fp: #here for debugging
    #     json.dump(detail_list, fp, indent=4)
    
    if wants_xlsx:

        try:
            xlsx_export(detail_list, cat_cfg["names"], output_path)
        except Exception as e: 
            print("Pls select a appropriate csv, if that does not help something else went wrong.") #just in case 


    plt.show()

            
def dic_initialiser(to_initialise):

    with open("./categories_configuration.json", 'r') as config: 
        cat_cfg = json.load(config)


    for elem in cat_cfg["names"]:
        to_initialise[elem] = 0


    return to_initialise


def group_lister(csv): #redundant a bit

    group_dic = {}

    for line in csv:
        
        len_csv = len(line) - 2

        missing_gen = 0
        missing_gen = line.count(' - , -')

        res = next(val for val in line if val != ' - , -') #how does it work?
        
        EDGAR_ID = res.split(",")[0]
        
            
        #core
        if len_csv - missing_gen == len_csv:
            group_dic[EDGAR_ID] = 0
        #singel
        elif len_csv - missing_gen == 1:
            group_dic[EDGAR_ID] = 1
        #disp
        else:
            group_dic[EDGAR_ID] = 2

    
    return group_dic


def dict_merger(EDGAR_ID, KEGG_ID, group_dic, database, singelton_dic, dispensable_dic, core_dic, KOFAM_file_path, detail_list, just_first): # questionable function, didnt get better with time, and even worse shield eyes
    
    #get edgar ids here
    
    try:
        if KEGG_ID != None:

            if just_first:

                if group_dic[EDGAR_ID] == 0:
                        if database[KEGG_ID]["categorie_sub"][0] in core_dic:
                            core_dic[database[KEGG_ID]["categorie_sub"][0]] += 1
                            save_details(KOFAM_file_path, 0, database[KEGG_ID]["categorie_sub"][0], database, KEGG_ID, detail_list, EDGAR_ID)
                        else:
                            core_dic['Other'] += 1
                            save_details(KOFAM_file_path, 0, 'Other', database, KEGG_ID, detail_list, EDGAR_ID)

                elif group_dic[EDGAR_ID] == 1:
                        if database[KEGG_ID]["categorie_sub"][0] in singelton_dic:
                            singelton_dic[database[KEGG_ID]["categorie_sub"][0]] +=1
                            save_details(KOFAM_file_path, 1, database[KEGG_ID]["categorie_sub"][0], database, KEGG_ID, detail_list, EDGAR_ID)
                        else:
                            singelton_dic['Other'] += 1
                            save_details(KOFAM_file_path, 1, 'Other', database, KEGG_ID, detail_list, EDGAR_ID)

                elif group_dic[EDGAR_ID] == 2:
                        if database[KEGG_ID]["categorie_sub"][0] in dispensable_dic:
                            dispensable_dic[database[KEGG_ID]["categorie_sub"][0]] += 1
                            save_details(KOFAM_file_path, 2, database[KEGG_ID]["categorie_sub"][0], database, KEGG_ID, detail_list, EDGAR_ID)
                        else:
                            dispensable_dic['Other'] += 1
                            save_details(KOFAM_file_path, 2, 'Other', database, KEGG_ID, detail_list, EDGAR_ID)
                            
            else:

                if group_dic[EDGAR_ID] == 0:
                    for elem in database[KEGG_ID]["categorie_sub"]:
                        if elem in core_dic:
                            core_dic[elem] += 1
                            save_details(KOFAM_file_path, 0, elem, database, KEGG_ID, detail_list, EDGAR_ID)
                        else:
                            core_dic['Other'] += 1
                            save_details(KOFAM_file_path, 0, 'Other', database, KEGG_ID, detail_list, EDGAR_ID)

                elif group_dic[EDGAR_ID] == 1:
                    for elem in database[KEGG_ID]["categorie_sub"]:
                        if elem in singelton_dic:
                            singelton_dic[elem] +=1
                            save_details(KOFAM_file_path, 1, elem, database, KEGG_ID, detail_list, EDGAR_ID)
                        else:
                            singelton_dic['Other'] += 1
                            save_details(KOFAM_file_path, 1, 'Other', database, KEGG_ID, detail_list, EDGAR_ID)

                elif group_dic[EDGAR_ID] == 2:
                    for elem in database[KEGG_ID]["categorie_sub"]:
                        if elem in dispensable_dic:
                            dispensable_dic[elem] += 1
                            save_details(KOFAM_file_path, 2, elem, database, KEGG_ID, detail_list, EDGAR_ID)
                        else:
                            dispensable_dic['Other'] += 1
                            save_details(KOFAM_file_path, 2, 'Other', database, KEGG_ID, detail_list, EDGAR_ID)

        else:

            if group_dic[EDGAR_ID] == 0:
                core_dic['Not Annotated'] += 1

            elif group_dic[EDGAR_ID] == 1:
                singelton_dic['Not Annotated'] +=1

            elif group_dic[EDGAR_ID] == 2:
                dispensable_dic['Not Annotated'] += 1

    except Exception as e: 
        print(f"{e} could not be found for some reason, check csv/txt, (or change DB)")



def save_details(KOFAM_file_path, pan_core_from, categorie, database, KEGG_ID, detail_list, EDGAR_ID): #uff
  
    if KOFAM_file_path not in detail_list:
        detail_list[KOFAM_file_path] = {}

    if pan_core_from not in detail_list[KOFAM_file_path]:
        detail_list[KOFAM_file_path][pan_core_from] = {}
    
    if categorie not in detail_list[KOFAM_file_path][pan_core_from]:
        detail_list[KOFAM_file_path][pan_core_from][categorie] = []

    detail_list[KOFAM_file_path][pan_core_from][categorie].append(f"{EDGAR_ID} ___ {database[KEGG_ID]['details']}")


def xlsx_export(detail_list, order, xlsxname):

    now = datetime.now()
    dt_string = now.strftime("%H_%M_%S__%d_%m")
    
    if xlsxname == None or xlsxname == "":
        xlsxname = dt_string

        for key in detail_list.keys():
            nameospath = os.path.basename(key)
            
            xlsxname+=f"_{nameospath[:3]}"

    workbook = xlsxwriter.Workbook(f'{xlsxname}.xlsx')
    
    row = 0
    page= 0
    
    cell_format = workbook.add_format({'bold': True})
    
    for fp in detail_list:

        page+=1
        row = 0
        nameoshere = os.path.basename(fp)
        worksheet = workbook.add_worksheet(f"{nameoshere[:3]}.({page})")

        worksheet.write(row, 0, fp, cell_format)
        row += 1
        worksheet.write(row, 0, "")
        row += 1

        for pc in detail_list[fp]: # translate 0123

            fraction=""

            if pc == 0:
                fraction = "Core Genes"
            elif pc == 1:
                fraction = "Singelton Genes"
            elif pc == 2:
                fraction = "Dispensable Genes"
            elif pc == 3:
                fraction = "Not Ordered into pan- dispensable- and core- genes"

            row += 2
            worksheet.write(row, 0, fraction, cell_format)
            row += 1
            worksheet.write(row, 0, "")
            row += 1

            for c in order:
                if c in detail_list[fp][pc]:
                    worksheet.write(row, 0, "")
                    row += 2
                    worksheet.write(row, 0, c, cell_format)
                    row += 1
                    worksheet.write(row, 0, "")
                    row += 1

                    for s in detail_list[fp][pc][c]:
                        worksheet.write(row, 0, s)
                        row += 1
    
    workbook.close()
    print(f"{xlsxname} | data exported")

def onclick_xlsx_export(order, num, filename, gui_on, detail_list, KOFAM_file_path, close): #just on click no confirmation needed, feels better, "gui on" add one! confirmation, order for v2 

    row = 0

    now = datetime.now()
    dt_string = now.strftime("%H_%M_%S__%d_%m")
    workbook = xlsxwriter.Workbook(f'{dt_string}_onclick.xlsx')
    cell_format = workbook.add_format({'bold': True})

    if KOFAM_file_path == None and order != 'Not Annotated':

        num = 3
        nameoshere = os.path.basename(filename)
        worksheet = workbook.add_worksheet(f"{nameoshere[:3]}.")

        worksheet.write(row, 0, filename, cell_format)
        row += 2
        worksheet.write(row, 0, "Not Ordered into pan- dispensable- and core- genes", cell_format)
        row += 2
        worksheet.write(row, 0, order, cell_format)
        row += 2
        worksheet.write_column(row, 0, detail_list[filename][num][order])

        workbook.close()
        print(f"{order}, {filename}, ({KOFAM_file_path}) | data exported")

    elif KOFAM_file_path == None and order == 'Not Annotated':
        print("It is not possible to export Gen-descriptions for Genes that were not annotated")

    if KOFAM_file_path != None and order != 'Not Annotated':

        if num == 0:
            fraction = "Core Genes"
        elif num == 2:
            fraction = "Singelton Genes"
        elif num == 1:
            fraction = "Dispensable Genes"

        switched = False

        if num == 1: #reversing the twist
            num = 2
            switched = True

        if num == 2 and not switched:
            num = 1
        
        nameoshere = os.path.basename(KOFAM_file_path)
        worksheet = workbook.add_worksheet(f"{nameoshere[:3]}.")

        worksheet.write(row, 0, KOFAM_file_path, cell_format)
        row += 2
        worksheet.write(row, 0, fraction, cell_format)
        row += 2
        worksheet.write(row, 0, order, cell_format)
        row += 2

        worksheet.write_column(row, 0, detail_list[KOFAM_file_path][num][order])

        workbook.close()
        print(f"{order}, {filename}, ({KOFAM_file_path}) | data exported")

    elif KOFAM_file_path != None and order == 'Not Annotated':
        print("It is not possible to export Gen-descriptions for Genes that were not annotated")


def remove_duplicates(input_path, output_path):

    set_non_unique_found = 0
    set_write = 0

    if output_path == None:
        output_path = "default_cut_faa_name.faa"
    try:
        unique = []

        p = re.compile(r'[>][\w]*[\s]{3}')

        with open(input_path, "r") as EDGAR_file:
            with open (output_path, "w") as nodupli_file:

                for line in EDGAR_file:
                    
                    if line.startswith('>'):

                        result = p.search(line)

                        if result == None: #fallback in case of ncbi faa, they also don't have duplicates either way (at least those i checked)
                            result = line

                        try: # looks horrific speed wise but should be fine for most files makes it more compatible with other faa files
                            identifier = (result.group(0))
                        except:
                            identifier = line

                        if identifier not in unique:
                            unique.append(identifier)
                            set_write = 1
                            
                        else:
                            print("Non unique found")
                            set_write = 0
                            set_non_unique_found = 1

                    if set_write == 1:
                        nodupli_file.write(line)

        if set_non_unique_found == 0: #hm botched
            os.remove(output_path)
            print("No duplicates found, no output file created")

        return set_non_unique_found

    except Exception as e: 
        print(e)
        print("Please select an (EDGAR) faa. file with duplicate IDs")


def auto_install_db(brite_type):

    if brite_type == None or brite_type == "":
        brite_type = "ko"

    print("Downloading and installing KEGG DB (takes 30 ish seconds)")
    resp = requests.get(f"https://www.kegg.jp/kegg-bin/download_htext?htext={brite_type}00001.keg&format=htext&filedir=")

    
    print("If this does not work check your brite type key or KEGG changed their download link / format")

    with open('rawdb.txt', 'wb') as f:
        f.write(resp.content)

    print("DB Installed, pls check DB validity")

    KEGG_to_json('rawdb.txt', "keggdb.json")
    update_dbpath("./keggdb.json")


def KEGG_to_json(path, output_path): #  there would have been https://www.kegg.jp/kegg-bin/download_htext?htext=ko00001.keg&format=json&filedir=, but its a little strange, not sure if useful

    try:
        
        f = open(path, "r")
        

        shift_for_map_a_to_b = ["Genetic Information Processing","Environmental Information Processing","Cellular Processes","Organismal Systems","Human Diseases"]
        #bad

        #regexp = re.compile(r'[D]\s{6}[K][0-9]{5}') #old one prevented buc
        regexp = re.compile(r'[K]{1}[0-9]{5}') # a littel more unspecific but works with all
        #over engineered?

        K_dic = {}
        A_for_dic = ""
        B_for_dic = ""

        for line in f:
            
            
            if line[:1] == "A":
                A_for_dic = line[7:-1]
                

            if line[:1] == "B":
                B_for_dic = line[9:-1]
                

            if A_for_dic in shift_for_map_a_to_b:
                B_for_dic = A_for_dic
                
                        
            if regexp.search(line): # line[:1] == "D": is as fast


                K_dic_key_value =	{
                                "categorie": [A_for_dic],
                                "categorie_sub": [B_for_dic],
                                "details": line[13:-1].strip() # 15:-1 on ref db, even with shoter ID this should include all detail 
                                
                }
                
                result = regexp.search(line)
                K_dic_key = result.group(0) # experimental! should work
             
                #K_dic_key = line[7:12] # [7:13] on ref db
                  
                if K_dic_key in K_dic:

                    if A_for_dic not in K_dic[K_dic_key]['categorie']:
                        K_dic[K_dic_key]['categorie'].append(A_for_dic)
                    if B_for_dic not in K_dic[K_dic_key]['categorie_sub']:
                        K_dic[K_dic_key]['categorie_sub'].append(B_for_dic)

                else:

                    K_dic[K_dic_key] = K_dic_key_value

        if output_path != None:
            with open(output_path, 'w') as fp:
                json.dump(K_dic, fp, indent=4)      
        else:
            with open(f'{path}.json', 'w') as fp:
                json.dump(K_dic, fp, indent=4)

    except Exception as e: 
        print(e)
        print("Pls select origninal KEGG Brite Hierarchy. Download: https://www.kegg.jp/kegg-bin/get_htext?ko00001.keg")


def remove_d_KO(path):

    K_number_list = []

    with open(path, "r") as EDGAR_file:
        with open (f"{path}_noduli_KO.txt", "w") as nodupli_file:
            for line in EDGAR_file:

                splitline = line.split("\t")

                if len(splitline) > 1:
                    K_number = splitline[1].strip("\n")

                if K_number not in K_number_list:
                    K_number_list.append(K_number)
                    nodupli_file.write(line)
                    
                else:
                    pass
                    #print("debug_message") #teststuff
                    
    print("Duplicate K numbers removed (and all \"not annotated\") Important Info: Also removes Gen \"duplicates\" within a single organism's set, if K number is deemed to be identical!")


def GUI():

    def explorer_win():

        app.setEntry("File 1", "")
        app.setEntry("File 2", "")
        app.setEntry("File 3", "")
        app.setEntry("Split CSV", "")
        box_entry = app.openBox(title=None, dirName=None, fileTypes=None, asFile=False, parent="Metabolomics_window", multiple=True, mode='r')

        for entry in box_entry: 
            if entry[-4:] == "json":
                app.setEntry("DB Path", str(entry))
                update_dbpath(entry)
                box_entry.remove(entry)

        for entry in box_entry: #need 2 seperate loops because of remove me thinks?
            if entry[-3:] == "csv":
                app.setEntry("Split CSV", str(entry))
                box_entry.remove(entry)
                
        try:
            app.setEntry("File 1", str(box_entry[0]))
            app.setEntry("File 2", str(box_entry[1]))
            app.setEntry("File 3", str(box_entry[2]))

        except:
            pass # 10/10 not all will be filled 


    def explorer_win_generic(calling_parent, to_fill):
        box_entry = app.openBox(title=None, dirName=None, fileTypes=None, asFile=False, parent="Metabolomics_window", multiple=False, mode='r')
        app.setEntry(to_fill, str(box_entry))


    def metabolomics_gui_analyse(button):

        wants_pancore = app.getCheckBox("split_box")

        if button == "Select Files":
            explorer_win()

        elif button == "Analyse" and wants_pancore == False: #was != True uh
            
            dbjson = app.getEntry("DB Path")
            file1 = app.getEntry("File 1")
            file2 = app.getEntry("File 2")
            file3 = app.getEntry("File 3")
            not_annotated = app.getCheckBox("not_annotated_box")
            pie = app.getCheckBox("pie_chart")
            wants_xlsx = app.getCheckBox("wants_xlsx")
            output_path = app.getEntry("Xlsx Name")
            app.hideSubWindow("onclick_print_window", useStopFunction=False)

            if not pie:
                app.showSubWindow("onclick_print_window")

            percent = app.getCheckBox("wants_percent")
            just_first = app.getCheckBox("just_first_box")

            analyse(dbjson, [file1, file2, file3], not_annotated, pie, output_path, wants_xlsx, percent, just_first, True, app)

        elif button == "Analyse" and wants_pancore == True:

            dbjson = app.getEntry("DB Path")
            file1 = app.getEntry("File 1")
            csv_path = app.getEntry("Split CSV")
            not_annotated = app.getCheckBox("not_annotated_box")
            pie = app.getCheckBox("pie_chart")
            wants_xlsx = app.getCheckBox("wants_xlsx")
            output_path = app.getEntry("Xlsx Name")
            app.hideSubWindow("onclick_print_window", useStopFunction=False)

            if not pie:
                app.showSubWindow("onclick_print_window")

            percent = app.getCheckBox("wants_percent")
            just_first = app.getCheckBox("just_first_box")
 
            pan_core_analysis_meta(csv_path, dbjson, file1, not_annotated, pie, output_path, wants_xlsx, percent, just_first, True, app)


    def genome_gui_analyse(button):

        if button == "Select CSV":
            explorer_win_generic("Pan-Core-Genome_window","CSV Path")
        elif button == "Analyse CSV":
            
            csvpath = app.getEntry("CSV Path")
            pie = app.getCheckBox("pie_csv")
            pan_core_analysis(csvpath, pie)


    def create_database_gui(button):
        
        if button == "Select KEGG DB":
            explorer_win_generic("Create Database_window", "KEGG DB Path")
        elif button == "Create json":
            
            KEGG_db_path = app.getEntry("KEGG DB Path")
            KEGG_db_path_name = app.getEntry("JSON DB Path")

            KEGG_to_json(KEGG_db_path, f"{KEGG_db_path_name}.json")


    def remove_duplicates_gui(button):
        
        if button == "Select EDGAR File":
            explorer_win_generic("Remove Duplicates_window", "EDGAR Path")
        elif button == "Rem. Duplicates":
            
            filtered_EDGAR_name = app.getEntry("Filtered EDGAR")
            EDGAR_path = app.getEntry("EDGAR Path")
            

            remove_duplicates(EDGAR_path, f"{filtered_EDGAR_name}.faa")

    
    def confirmation_window_gui(button):

        if button == "Yes":
            #app.thread(auto_install_db) #Not wroth the risk, it's just like 15 sek of freeze anyway, app jar multithreading unstable?
            Brite_type = app.getEntry("Brite Type")
            auto_install_db(Brite_type)

            app.infoBox("Installed","DB install process finished, \n please check DB validity", parent="Confirmation_window") #yes this not a real check

               
        elif button == "No":
            app.hideSubWindow("Confirmation_window", useStopFunction=False)


    def press_main_menu(button):

        if button == "Metabolic-Profiling":
            print("KofamKoala outputs for File 1,2,3 | for direct processing of .faa's see -h")
            dbjson = load_dbpath()
            app.setEntry("DB Path", str(dbjson))
            app.showSubWindow("Metabolomics_window")
        elif button == "Pan-Core-Genome":
            app.showSubWindow("Pan-Core-Genome_window")
        elif button == "Create Database":
            app.showSubWindow("Create Database_window")
        elif button == "Remove Duplicates":
            app.showSubWindow("Remove Duplicates_window")
        elif button == "Install DB":
            app.showSubWindow("Confirmation_window")

    app = gui("CoGeMPA")
    app.setFont(14)
    app.setSize(800,400)
    app.setBg("gray")
    app.addLabel("title", "Main Menu")
    app.getLabelWidget("title").config(font="Arial 30 bold underline")

    app.addButtons(["Metabolic-Profiling" ,"Pan-Core-Genome", "Create Database", "Remove Duplicates", "Install DB"], press_main_menu)


    def checkStop():
        app.hideAllSubWindows(useStopFunction=False)
        plt.close('all')
        return app.yesNoBox("Confirm Exit", "Are you sure you want to exit the application?")


    app.setStopFunction(checkStop)
    
    #sub window sure
    app.startSubWindow("Confirmation_window", title="Confirmation")
    app.setSize(600,400)
    app.setBg("gray")
    app.addLabel("text1", "Are you sure you want to download and install a new database?")
    app.addLabelEntry("Brite Type")
    app.addLabel("text2", "(This takes about 30 seconds, the GUI is unresponsive for this time.)")
    app.addLabel("brite_text", "(Brite type uses \"Reference\" if not set.)")
    app.addWebLink("Brite Type Names (in blue)", "https://www.kegg.jp/kegg-bin/list_menu?id=ko00001.keg")
    app.addButtons(["Yes","No"], confirmation_window_gui)

    app.config_window_shown = False
    #app.setOnTop(stay=True) annoying
    app.stopSubWindow()
    
    #sub window Metabolomics
    app.startSubWindow("Metabolomics_window", title="Metabolic Profiling") # 'Metabolomics' misnomer
    app.setSize(600,400)
    app.setBg("gray")

    app.addLabelEntry("DB Path")

    app.addLabelEntry("File 1")
    app.addLabelEntry("File 2")
    app.addLabelEntry("File 3")
    app.addLabelEntry("Split CSV")
    app.addLabelEntry("Xlsx Name")
    app.addNamedCheckBox("Show not Annotated","not_annotated_box")
    app.addNamedCheckBox("Pie Chart","pie_chart")
    app.addNamedCheckBox("Just show 1'st category","just_first_box")
    app.addNamedCheckBox("Show \"PanCore\" Distribution","split_box")
    app.addNamedCheckBox("Export Data","wants_xlsx")
    app.addNamedCheckBox("Show Barchart in %","wants_percent")
    app.addButtons(["Select Files","Analyse"], metabolomics_gui_analyse)
    app.config_window_shown = False
    #app.setOnTop(stay=True)
    app.stopSubWindow()

    #sub window Pan-Core-Genome
    app.startSubWindow("Pan-Core-Genome_window", title="Pan-Core-Genome")
    app.setSize(600,400)
    app.setBg("gray")
    app.addLabelEntry("CSV Path")
    app.addNamedCheckBox("Pie Chart","pie_csv")
    app.addButtons(["Select CSV","Analyse CSV"], genome_gui_analyse)

    app.config_window_shown = False
    #app.setOnTop(stay=True)
    app.stopSubWindow()

    #sub window Create Database
    app.startSubWindow("Create Database_window", title="Create Database")
    app.setSize(600,400)
    app.setBg("gray")
    app.addLabelEntry("KEGG DB Path")
    app.addLabelEntry("JSON DB Path")
    app.addButtons(["Select KEGG DB","Create json"], create_database_gui)
    app.config_window_shown = False
    #app.setOnTop(stay=True)
    app.stopSubWindow()

    #sub window Remove Duplicates
    app.startSubWindow("Remove Duplicates_window", title="Remove Duplicates")
    app.setSize(600,400)
    app.setBg("gray")
    app.addLabelEntry("EDGAR Path") #had name at end 
    app.addLabelEntry("Filtered EDGAR") #had name at end
    app.addButtons(["Select EDGAR File","Rem. Duplicates"], remove_duplicates_gui)
    app.addLabel("infolabel_2", "Optional any faa file can be processed with this")
    app.config_window_shown = False
    #app.setOnTop(stay=True)
    app.stopSubWindow()

    #sub window onclick_print
    app.startSubWindow("onclick_print_window", title="Info")
    app.setSize(600,200)
    app.addLabel("infolabel", "To export data of a specific graph click it")
    app.config_window_shown = False
    #app.setOnTop(stay=True)
    app.stopSubWindow()

    app.go()


def update_dbpath(newpath):

    with open("./categories_configuration.json", 'r') as config: 
        cfg = json.load(config)

    cfg["dbpath"] = newpath
    
    with open("./categories_configuration.json", 'w') as config: 
        json.dump(cfg, config, indent=4)


def load_dbpath():
    with open("./categories_configuration.json", 'r') as config: 
                cfg = json.load(config)
    return(cfg["dbpath"])


def KOFAM_subprocess(args):

    to_process= [args.file1, args.file2, args.file3]
    to_process = list(filter(None, to_process)) 

    processed= []

    if args.rem_auto:
        for i, e in enumerate(to_process):

            set_non_unique_found = remove_duplicates(e, f"{e}_noduplicates.faa")
            
            if set_non_unique_found == 1:
                to_process[i] = f"{e}_noduplicates.faa"
                print("Duplicates in EDGAR faa. removed")
            elif set_non_unique_found == 0:
                to_process[i] = e

    print(to_process)
    print("Files that will be processed need to be in the .faa format (except already preprocessed Kofam Koala outputs)")
    print("This can take quite long")

    for f in to_process:
        if f != None and f[-4:] == ".faa":
            print(f"Now processing {f}")
            p1 = subprocess.Popen(["ruby", "exec_annotation" , "-f", "mapper", "-o", f"{f}_KOFAM_output.txt", f]) 
            processed.append(f"{f}_KOFAM_output.txt")                       
            p1.wait()
            print(f"One faa. file finished ({f})")

        else:
            processed.append(f)
            print(f"One file was already preprocessed ({f})")

    if len(to_process) == 2:
        processed.append(None)

    elif len(to_process) == 1:
        processed.append(None)
        processed.append(None)

    print(processed)
    args.file1, args.file2, args.file3 = processed[0], processed[1], processed[2]


def main():
    args = argp()

    if args.fullauto:
        KOFAM_subprocess(args)
    
    if args.start_gui:
        GUI()

    elif args.met_analysis:
        
        try:
            if args.dbjson == None:
                dbjson = load_dbpath()
                analyse(dbjson, [args.file1, args.file2, args.file3], args.not_annotated, args.pie_chart, args.output, args.export_excel, args.percent, args.just_first, False, None)
        

            else:
                update_dbpath(args.dbjson)
                analyse(args.dbjson, [args.file1, args.file2, args.file3], args.not_annotated, args.pie_chart, args.output, args.export_excel, args.percent, args.just_first, False, None)

        except:
            print("Please Select an appropriate Database Path with -dbj")

    elif args.csv_analysis:
        pan_core_analysis(args.csvpath, args.pie_chart)

    elif args.rem_duplicates:
        remove_duplicates(args.EDGAR, args.output)

    elif args.make_database:
        KEGG_to_json(args.database, args.output)


    elif args.csv_meta:
        try:
            if args.dbjson == None:
                dbjson = load_dbpath()
                pan_core_analysis_meta(args.csvpath, dbjson, args.file1, args.not_annotated, args.pie_chart, args.output, args.export_excel, args.percent, args.just_first, False, None)

            else:
                update_dbpath(args.dbjson)
                pan_core_analysis_meta(args.csvpath, args.dbjson, args.file1, args.not_annotated, args.pie_chart, args.output, args.export_excel, args.percent, args.just_first, False, None)
        except Exception as e: 
            print(e)
            print("Please Select appropriate files")

    elif args.install:
        auto_install_db(args.database)

    elif args.remko:
        remove_d_KO(args.file1)

    else:
        print("Please select a function, try -h for help")


    
if __name__ == "__main__":
    main()