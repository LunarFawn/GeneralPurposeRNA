import mysql.connector
from mysql.connector import Error
import pandas as pd

import nupackAPI_Sara2_Ver1 as nupackAPI
import Sara2_API_Python3 as sara2

def Create_server_connectionCursor(host_name, user_name, user_password):
    connection = None
    cursor = None
    try:
        connection = mysql.connector.connect(
            host=host_name,
            user=user_name,
            passwd=user_password
        )
        cursor = connection.cursor()
        print("MySQL Database connection successful")
    except Error as err:
        print(f"Error: '{err}'")

    return connection, cursor 

def Execute_query(cursor, connection, query):
    isSuccess = False
    try:
        cursor.execute(query)
        connection.commit()
        print("Query successful")
        isSuccess=True
    except Error as err:
        print(f"Error: '{err}'")
    return isSuccess

def MakeDB(cursor, databaseName):
    query = "CREATE DATABASE {database};".format(database=databaseName)
    isMade = False
    try:
        cursor.execute(query)
        print("Database created successfully")
        isMade=True
    except Error as err:
        print(f"Error: '{err}'")
    return isMade

def MakeTable(cursor, connection, tableName):
    #make a empty table with only a ID column that also acts as teh primary key.
    #will add collumns as needed
    query = "CREATE TABLE {tableName} (id INT AUTO_INCREMENT PRIMARY KEY);".format(tableName=tableName)
    isMade = Execute_query(cursor, connection, query)
    return isMade

def InsertCollumn(cursor, connection, tableName, columnName, datatype):
    query = "ALTER TABLE {tableName} ADD {columnName} {datatype};".format(tableName=tableName, columnName=columnName)

def GenerateSQL(puzzle: sara2.puzzleData, databaseName):   
    connection, cursor = Create_server_connectionCursor("localhost", "rnauser", "rna")
    dbMade = MakeDB(cursor, databaseName)

    if dbMade is True:
        #now make a table for each part of each desing namedtuple list
        #there will also be a table for the design info and then each nucleotide table will be an embeded table
        #that will have each nuc pair possible with a value for prob if present in fold or null/something if not        
        #make the pairs table first
        pairsTable = "PAIRS"
        tableMade = MakeTable(cursor, connection, pairsTable)
        
        if tableMade is True:
            pass

