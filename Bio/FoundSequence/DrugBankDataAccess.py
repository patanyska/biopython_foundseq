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

    Variables: - idication: is the patology or problem founded sequence
"""
def get_drug_by_indication(connection,indication):
   try:
       strings_indication=[]
       splited_indication=indication.split()
       composed_indication=""
       pos=0
       for i in splited_indication:
           if(pos==0):
               strings_indication.append(i)
               composed_indication=i
           else:
               strings_indication.append(composed_indication+" "+i)
               composed_indication=composed_indication+" "+i
           pos=pos+1
               
       indication=""
       counti=1
       for s in strings_indication:
           if(counti==1):
               indication=indication +"(d.indication like '%"+s+"%' OR "
           else:
               if(counti<len(strings_indication)):
                   indication=indication +"d.indication like '%"+s+"%' OR "
               else:
                   indication=indication +"d.indication like '%"+s+"%') COLLATE NOCASE"
           counti=counti+1
       
       pathway=""
       counti=1
       for s in strings_indication:
           if(counti==1):
               pathway=pathway +"(pt.name like '%"+s+"%' OR "
           else:
               if(counti<len(strings_indication)):
                   pathway=pathway +"pt.name like '%"+s+"%' OR "
               else:
                   pathway=pathway +"pt.name like '%"+s+"%') COLLATE NOCASE"
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
                       ("""+indication+""" OR """+pathway+""") 
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

