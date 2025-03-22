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
def get_drug_by_diasease(connection,disease):
    records=[]
    try:

        words_disease=[]
        splited_words=disease.split(" ") #split each disease in words ex:[heart,attack]
        composed_disease=""
        pos=0
        for i in splited_words:
            if(pos==0):
                words_disease.append(i)
                composed_disease=i
            else:
                words_disease.append(composed_disease+" "+i)
                composed_disease=composed_disease+" "+i
            pos=pos+1

            
        if(len(words_disease)==1):
            col_description="(d.description like '%"+words_disease[0]+"%') COLLATE NOCASE"
            col_indication="(d.indication like '%"+words_disease[0]+"%') COLLATE NOCASE"
            col_pathway="(pt.name like '%"+words_disease[0]+"%') COLLATE NOCASE"
        else:
            col_description=""
            counti=1
            for word in words_disease:
                if(counti==1):
                    col_description=col_description +"(d.description like '%"+word+"%' OR "
                else:
                    if(counti<len(words_disease)):
                        col_description=col_description +"d.description like '%"+word+"%' OR "
                    else:
                        col_description=col_description +"d.description like '%"+word+"%') COLLATE NOCASE"
                counti=counti+1
                    
            col_indication=""
            counti=1
            for word in words_disease:
                if(counti==1):
                    col_indication=col_indication +"(d.indication like '%"+word+"%' OR "
                else:
                    if(counti<len(words_disease)):
                        col_indication=col_indication +"d.indication like '%"+word+"%' OR "
                    else:
                        col_indication=col_indication +"d.indication like '%"+word+"%') COLLATE NOCASE"
                counti=counti+1
            
            col_pathway=""
            counti=1
            for word in words_disease:
                if(counti==1):
                    col_pathway=col_pathway +"(pt.name like '%"+word+"%' OR "
                else:
                    if(counti<len(words_disease)):
                        col_pathway=col_pathway +"pt.name like '%"+word+"%' OR "
                    else:
                        col_pathway=col_pathway +"pt.name like '%"+word+"%') COLLATE NOCASE"
                counti=counti+1
            
        cursor=connection.cursor()
        sql="""SELECT distinct(d.name)            
        ,d.drugbank_id
        ,d.description
        ,d.state
        ,d.indication
        ,GROUP_CONCAT(distinct(p.route)) as route                          
        ,GROUP_CONCAT(distinct(p.country)) as country                    
        FROM drug d 
            INNER JOIN groups g on d.drugbank_id=g.drugbank_id                   
            INNER JOIN products p on d.drugbank_id=p.drugbank_id
            INNER JOIN pathaway pt on d.drugbank_id=pt.drugbank_id   
                WHERE 
                    ("""+col_indication+""") 
                    AND p.approved='true'
                    AND (p.ended_marketing_on IS "" OR p.ended_marketing_on>=date('now'))
                    AND d.drugbank_id NOT IN (select drugbank_id from groups gr where gr.description like "withdrawn")
                GROUP BY d.drugbank_id
                HAVING d.name IS NOT NULL
        UNION ALL

        SELECT distinct(d.name)            
        ,d.drugbank_id
        ,d.description
        ,d.state
        ,d.indication
        ,GROUP_CONCAT(distinct(p.route)) as route                          
        ,GROUP_CONCAT(distinct(p.country)) as country                    
        FROM drug d 
            INNER JOIN groups g on d.drugbank_id=g.drugbank_id                   
            INNER JOIN products p on d.drugbank_id=p.drugbank_id
            INNER JOIN pathaway pt on d.drugbank_id=pt.drugbank_id   
                WHERE 
                    ("""+col_pathway+""") 
                    AND p.approved='true'
                    AND (p.ended_marketing_on IS "" OR p.ended_marketing_on>=date('now'))
                    AND d.drugbank_id NOT IN (select drugbank_id from groups gr where gr.description like "withdrawn")
                GROUP BY d.drugbank_id
                HAVING d.name IS NOT NULL
        UNION ALL

        SELECT distinct(d.name)            
        ,d.drugbank_id
        ,d.description
        ,d.state
        ,d.indication
        ,GROUP_CONCAT(distinct(p.route)) as route                          
        ,GROUP_CONCAT(distinct(p.country)) as country                    
        FROM drug d 
            INNER JOIN groups g on d.drugbank_id=g.drugbank_id                   
            INNER JOIN products p on d.drugbank_id=p.drugbank_id
            INNER JOIN pathaway pt on d.drugbank_id=pt.drugbank_id   
                WHERE 
                    ("""+col_description+""") 
                    AND p.approved='true'
                    AND (p.ended_marketing_on IS "" OR p.ended_marketing_on>=date('now'))
                    AND d.drugbank_id NOT IN (select drugbank_id from groups gr where gr.description like "withdrawn")
                GROUP BY d.drugbank_id
                HAVING d.name IS NOT NULL
        """
        
        cursor.execute(sql)
        columns = [desc[0] for desc in cursor.description]  # Get column names
        records.extend([dict(zip(columns, row)) for row in cursor.fetchall()])  # Convert to dict list

        return records
    except:
       raise ValueError(
            f"A problem occurred during geting data from DrugBank")
  
   
#endregion

