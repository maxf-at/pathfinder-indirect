#!/usr/bin/env python3
# coding: utf-8

import sys
import random
import string
import RNA
import numpy as np
import os
import subprocess
from dataclasses import dataclass, field

import matplotlib.pyplot as plt
from PIL import Image
from IPython.display import SVG, display
from helper import print_moves

from sqlalchemy import create_engine
from sqlalchemy import Table, Column, Integer, String, MetaData, ForeignKey



class RNA_db:
    # http://docs.sqlalchemy.org/en/latest/core/engines.html
    DB_ENGINE = {
        "sqlite": 'sqlite:///{DB}'
    }

    # Main DB Connection Ref Obj
    db_engine = None
    def __init__(self, username='', password='', dbname=''):
        dbtype = "sqlite"
        dbtype = dbtype.lower()
        if dbtype in self.DB_ENGINE.keys():
            engine_url = self.DB_ENGINE[dbtype].format(DB=dbname)
            self.db_engine = create_engine(engine_url)


    def print_all_data(self, table='', query=''):
        query = query if query != '' else "SELECT * FROM '{}';".format(table)
        print("query:", query)
        with self.db_engine.connect() as connection:
            try:
                result = connection.execute(query)
            except Exception as e:
                print(e)
            else:
                for row in result:
                    print(row) # print(row[0], row[1], row[2])
                result.close()
        print("\n")



    def get_id(self, sequence, s1, s2, dataset):

        query = f"SELECT sequence, s1, s2 FROM sequences_with_structures\
                  WHERE sequence='{sequence}' AND s1='{s1}' AND s2='{s2}'"
        # self.print_all_data(query=query)

        bp_dist = RNA.bp_distance(s1, s2)
        length = len(sequence)

        result = []


        with self.db_engine.connect() as connection:
            try:
                result = connection.execute(query)
            except Exception as e:
                print(e)            
            data = result.fetchall()

            if not data:
                query = f"INSERT INTO sequences_with_structures (sequence, s1, s2, length, bp_dist, dataset)\
                          VALUES ('{sequence}','{s1}', '{s2}', '{length}', '{bp_dist}', '{dataset}');"
                try:
                    result = connection.execute(query)
                except Exception as e:
                    print(e)

            # get id
            query = f"SELECT id FROM sequences_with_structures\
                      WHERE sequence='{sequence}' AND s1='{s1}' AND s2='{s2}'"

            data = connection.execute(query)       
            result = data.fetchall()            
            id = result[0][0]

            return id
        
    def insert_path(self, sequence, s1, s2, dataset, moves, info):
        id = self.get_id(sequence, s1, s2, dataset)

        query = f"INSERT OR REPLACE INTO indirect_paths (moves, id, description)\
                          VALUES ('{moves}','{id}', '{info}');"

        # query = f"IF EXISTS(select * from indirect_paths where id={id} AND description='{info}')\
        #         update indirect_paths set moves='{moves}' where id={id}\
        #         ELSE\
        #         insert into test(name) values('john');"
        with self.db_engine.connect() as connection:
            try:
                connection.execute(query)
            except Exception as e:
                print(e)            



    def get_paths(self, dataset, description):

        query = f"SELECT sequence, s1, s2, sequences_with_structures.id, moves FROM sequences_with_structures\
                LEFT JOIN indirect_paths on indirect_paths.id = sequences_with_structures.id\
                WHERE dataset='{dataset}' AND description='{description}'\
                ORDER BY sequences_with_structures.id"
        # self.print_all_data(query=query)

        # bp_dist = RNA.bp_distance(s1, s2)
        # length = len(sequence)
        # result = []


        with self.db_engine.connect() as connection:
            try:
                result = connection.execute(query)
            except Exception as e:
                print(e)            
            data = result.fetchall()

            for i, (sequence, s1, s2, id, moves) in enumerate(data):
                moves = self.process_moves(moves)
                yield i, sequence, s1, s2, id, moves
                # print (sequence, s1, s2, id, moves)
            
            # print (data)

            # result = data.fetchall()

    def process_moves(self, move_strings):
        """
        input: "[(0, 0), (9, 21), (-35, -56)]" as string from DB
        output: list of tuples
        """

        # remove []
        move_strings = move_strings[1:-1]
        move_strings = move_strings.split("), (")
        
        moves = []
        for move in move_strings:
            i, j = move.replace('(', '').replace(')', '').split(",")
            i = int(i)
            j = int(j)
            moves.append((i,j))
        return moves


        # Sample Query Joining
        # query = "SELECT u.last_name as last_name, " \
        #         "a.email as email, a.address as address " \
        #         "FROM {TBL_USR} AS u " \
        #         "LEFT JOIN {TBL_ADDR} as a " \
        #         "WHERE u.id=a.user_id AND u.last_name LIKE 'M%';" \
        #     .format(TBL_USR=USERS, TBL_ADDR=ADDRESSES)
        # self.print_all_data(query=query)

if __name__ == '__main__':

    sequence = "UCACUGAGGCUUGUUCGCAAAUCACUGCAAUUAGAUAUGACUCACGAUAUGGGGCACGGUGCAUACAUAC"
    s1       = ".(((((.(.((((((((....((...............))....))).))))).).)))))........."
    s2       = ".(((((..((((..(((..(.(((.((........))))).)..)))....)))).)))))........."
    
    rna_db = RNA_db(dbname='./sample_db/rna_samples.sqlite')
     
    # insert example
    # moves = [(0,0),(1,2)]
    # info = "TABU"
    # database = "60_indirect"
    # rna_db.insert_path(sequence, s1, s2, database, moves, info)

    # dataset = "indirect_60.csv"
    # dataset = "indirect_ea_60.csv"
    dataset = "indirect_input_80.csv"
    # dataset = "indirect_ea_100.csv"

    description = "RNAeapath"
    # description = "TABU"
    
    for i, sequence, s1, s2, id, moves in rna_db.get_paths(dataset, description):
        
        print_moves(sequence, s1, s2, moves)
        
        # print (moves)



    # Create Tables
    # dbms.create_db_tables()

    # dbms.insert_single_data()
    # dbms.print_all_data(mydatabase.USERS)
    # dbms.print_all_data(mydatabase.ADDRESSES)
    # dbms.sample_query() # simple query
    # dbms.sample_delete() # delete data
    # dbms.sample_insert() # insert data
    # dbms.sample_update() # update data