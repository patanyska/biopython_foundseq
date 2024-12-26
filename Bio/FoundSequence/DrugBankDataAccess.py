# Copyright 2023 by Patricia Nogueira  All rights reserved.
#
# This file is part of my masters degree project
import os
import sqlite3
import datetime
import xml.etree.ElementTree as et
from pathlib import PurePath


    
def openConnection():
    mainpath=PurePath(__file__).parent
    db_path=os.path.join(str(mainpath) + "\\db_\\db_drugbank.db")
    connection=sqlite3.connect(db_path)
    return connection

def closeConnection(connection):
    connection.close()

#region Queries
"""
    Query all approved drugs for indication 

    Variables: - indication: is the patology or problem founded sequence
"""
def get_drug_by_diasease(connection,diasease):
    try:
        strings_disease=[]
        splited_disease=diasease.split()
        composed_disease=""
        pos=0
        for i in splited_disease:
           if(pos==0):
               strings_disease.append(i)
               composed_disease=i
           else:
               strings_disease.append(composed_disease+" "+i)
               composed_disease=composed_disease+" "+i
           pos=pos+1

        
        if(len(strings_disease)==1):
           col_description="(d.description like '%"+strings_disease[0]+"%') COLLATE NOCASE"
           col_indication="(d.indication like '%"+strings_disease[0]+"%') COLLATE NOCASE"
           col_pathway="(pt.name like '%"+strings_disease[0]+"%') COLLATE NOCASE"
        else:
            col_description=""
            counti=1
            for s in strings_disease:
                if(counti==1):
                    col_description=col_description +"(d.description like '%"+s+"%' OR "
                else:
                    if(counti<len(strings_disease)):
                        col_description=col_description +"d.description like '%"+s+"%' OR "
                    else:
                        col_description=col_description +"d.description like '%"+s+"%') COLLATE NOCASE"
                counti=counti+1
                    
            col_indication=""
            counti=1
            for s in strings_disease:
                if(counti==1):
                    col_indication=col_indication +"(d.indication like '%"+s+"%' OR "
                else:
                    if(counti<len(strings_disease)):
                        col_indication=col_indication +"d.indication like '%"+s+"%' OR "
                    else:
                        col_indication=col_indication +"d.indication like '%"+s+"%') COLLATE NOCASE"
                counti=counti+1
            
            col_pathway=""
            counti=1
            for s in strings_disease:
                if(counti==1):
                    col_pathway=col_pathway +"(pt.name like '%"+s+"%' OR "
                else:
                    if(counti<len(strings_disease)):
                        col_pathway=col_pathway +"pt.name like '%"+s+"%' OR "
                    else:
                        col_pathway=col_pathway +"pt.name like '%"+s+"%') COLLATE NOCASE"
                counti=counti+1
           
        cursor=connection.cursor()
        sql="""SELECT distinct(d.name)            
           ,d.drugbank_id
           ,d.description
           ,d.cas_number
           ,d.unii
           ,d.state
           ,d.indication
           ,d.pharmacodynamics
           ,d.mechanism_of_action          
           ,GROUP_CONCAT(distinct(p.name)) as product_name            
           ,GROUP_CONCAT(distinct(p.labeller)) as labeller                         
           ,GROUP_CONCAT(distinct(p.route)) as route                          
           ,GROUP_CONCAT(distinct(p.country)) as country                    
           FROM drug d 
               LEFT JOIN groups g on d.drugbank_id=g.drugbank_id                   
               LEFT JOIN products p on d.drugbank_id=p.drugbank_id
               LEFT JOIN pathaway pt on d.drugbank_id=pt.drugbank_id   
                   WHERE 
                       ("""+col_indication+""" OR """+col_pathway+""" OR """+col_description+""") 
                       AND p.approved='true'
                       AND (p.ended_marketing_on IS "" OR p.ended_marketing_on>=date('now'))
                       AND d.drugbank_id NOT IN (select drugbank_id from groups gr where gr.description like "withdrawn")"""
       
       
        cursor.execute(sql)
        records = cursor.fetchall()
        return records
    except:
       raise ValueError(
            f"A problem occurred during geting data from DrugBank")
  
   
#endregion

